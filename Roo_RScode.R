# 22 April 2025
# Run reproductive success model

## Set up ----------------------------------------------------------------------

# load packages
library(tidyverse)
library(lubridate)
library(beepr)
library(coda)
library(nimble)
library(foreach)
library(parallel)
library(doParallel)
registerDoParallel(3)

# load data
source("wrangleData_en.R")
enData <- wrangleData_en(dens.data = "data/WPNP_Methods_Results_January2025.xlsx",
                         veg.data  = "data/biomass data April 2009 - Jan 2025_updated Feb2025.xlsx",
                         wea.data  = "data/Prom_Weather_2008-2023_updated Jan2025 RB.xlsx",
                         wind.data = "data/POWER_Point_Daily_20080101_20241231_10M.csv")

source("wrangleData_rs.R")
rsData <- wrangleData_rs(rs.data = "data/RSmainRB_Mar25.xlsx",
                         obs.data = "data/PromObs_2008-2019.xlsx",
                         known.age = TRUE, cum.surv = TRUE, surv.sep1 = TRUE)

# check that all years are represented
setequal(1:17, unique(rsData$year)) # should be TRUE

# create Nimble lists
myData <-  list(rs   = rsData$survS1,
                id   = rsData$id,
                year = rsData$year,
                age  = rsData$age-2, # SO THAT AGE STARTS AT 1
                dens = enData$dens,
                veg  = enData$veg,
                win  = enData$win)

myConst <- list(n      = rsData$n,
                nIDr = rsData$nIDr,
                nYear  = rsData$nYear,
                nAge   = rsData$nAge,
                nAgeC  = rsData$nAgeC)
# TODO: deal with missing environment

# Switches/toggles
testRun <- TRUE # or FALSE


## Model -----------------------------------------------------------------------

myCode = nimbleCode({
  
  #### Likelihood & constraints ####
  # individual rs function
  for(x in 1:n){
    rs[x] ~ dbern(rsI[x])
    logit(rsI[x]) <- logit(Mu.rsI[age[x]]) +
      # BetaD.rs * dens[year[x]] +
      # BetaV.rs * veg[year[x]] +
      # BetaW.rs * win[year[x]] +
      EpsilonI.rsI[id[x]] +
      EpsilonT.rsI[year[x]]
  }
  
  # age-specific rs function
  # use parameters estimated from individual data above
  # to predict age-specific reproductive success (rsA) here!
  for(a in 1:nAge){
    for(t in 1:nYear){
      logit(rsA[a, t]) <- logit(Mu.rsA[a]) + # rsA becomes svPY!
        # BetaD.rs * dens[t] +
        # BetaV.rs * veg[t] +
        # BetaW.rs * win[t] +
        EpsilonT.rsA[t]
    }
  }
    
  ##### Priors ####
  # priors for fixed effects
  for(a in 1:nAge){
    Mu.rsI[a] ~ dunif(0, 1)
    Mu.rsA[a] ~ dunif(0, 1)
  }
  
  # Beta.dens ~ dunif(-2, 2) # could be dunif(-5, 5) if need be
  # Beta.veg  ~ dunif(-2, 2) # could be dunif(-5, 5) if need be
  # Beta.win  ~ dunif(-2, 2) # could be dunif(-5, 5) if need be
  
  # priors for random effects
  for(i in 1:nIDr){
    EpsilonI.rsI[i] ~ dnorm(0, sd = SigmaI.rsI)
  }
  
  for(t in 1:nYear){
    EpsilonT.rsI[t] ~ dnorm(0, sd = SigmaT.rsI)
    EpsilonT.rsA[t] ~ dnorm(0, sd = SigmaT.rsA)
  }
  
  # priors for sigma
  SigmaI.rsI ~ dunif(0, 100)
  SigmaT.rsI ~ dunif(0, 100)
  SigmaT.rsA ~ dunif(0, 100)
  
})


## Assemble --------------------------------------------------------------------

# to serialize
# create Nimble function
paraNimble <- function(mySeed, myCode, myConst, myData,
                       rsData = rsData, enData = enData, testRun){
  
  library(nimble)
  set.seed(mySeed)
  
  n = myConst$n
  nIDr = myConst$nIDr
  nYear = myConst$nYear
  nAge = myConst$nAge
  
  # assign initial values
  source("simulateInits.R")
  myInits <- simulateInits(n = n, nIDs = 0, nIDr = nIDr, nYear = nYear, nAge = nAge,
                           age = myData$age, dens = myData$dens, veg = myData$veg, win = myData$win)
                           # nNoVeg = myConst$nNoVeg, nNoWin = myConst$nNoWin
  
  # assemble model
  myMod <- nimbleModel(code = myCode,
                       data = myData,
                       constants = myConst,
                       inits = myInits
                       )
  
  # select parameters to monitor
  params = c(# RS model
             "Mu.rsI", "Mu.rsA",
             # "Beta.dens", "Beta.veg", "Beta.win",
             "EpsilonI.rsI", "EpsilonT.rsI", "EpsilonT.rsA",
             "SigmaI.rsI", "SigmaT.rsI", "SigmaT.rsA"
  )
  
  # MCMC settings
  if(testRun){
    niter   <- 10
    nburnin <- 0
    nthin   <- 1
  }else{
    niter   <- 10000
    nburnin <- 6000
    nthin   <- 1
  }
  
  cModel <- compileNimble(myMod)
  myMCMC <- buildMCMC(cModel, monitors = params)
  cmyMCMC <- compileNimble(myMCMC, project = myMod)
  samples <- runMCMC(cmyMCMC,
                     niter = niter,
                     nburnin = nburnin,
                     thin = nthin,
                     setSeed = mySeed,
                     samplesAsCodaMCMC = T,
                     summary = T)
  
  return(samples)
}


## Run model -------------------------------------------------------------------

# serialized
mySeed <- 1:3
nchains <- 3

start.t <- Sys.time()
this_cluster <- makeCluster(3)
samples <- parLapply(X = 1:3,
                     cl = this_cluster,
                     fun = paraNimble,
                     myCode = myCode,
                     myConst = myConst,
                     myData = myData,
                     testRun = testRun)

beep(sound = 2)
stopCluster(this_cluster)
dur = now() - start.t; dur

library(MCMCvis)
out.mcmc <- samples %>% map(~as.mcmc(.$samples)) %>% as.mcmc.list()

MCMCsummary(out.mcmc,
            params = c('Mu.rsI', 'Mu.rsA', 'SigmaI.rsI', 'SigmaT.rsI', 'SigmaT.rsA'),
            n.eff = TRUE, round = 2)

