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
enData <- wrangleData_env(dens.data = "data/WPNP_Methods_Results_January2025.xlsx",
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
                nID.rs = rsData$nID.rs,
                nYear  = rsData$nYear,
                nAge   = rsData$nAge,
                nAgeC  = rsData$nAgeC)
# TODO: deal with missing environment

# Switches/toggles
testRun <- TRUE # or FALSE


## Parameters ------------------------------------------------------------------

# n = number of observations, or reproductive events
# nID.rs = number of unique kangaroos in the dataset
# nYear = number of years in the dataset
# nAge = number of ages in the analysis (3 through 19 years old, so 17 ages)
# nAgeC = number of age classes in the analysis (not used so far in RS analysis)


## Model -----------------------------------------------------------------------

myCode = nimbleCode({
  
  #### First attempt ####
  # likelihood
  for (x in 1:n){
    rs[x] ~ dbern(rsI[id[x], year[x]])
  }

  # constraints
  for (i in 1:nID.rs){
    for (t in 1:nYear){
      logit(rsI[i, t]) <- logit(Mu.rsI[age[i, t]]) +
        # BetaD.rs * dens[t] +
        # BetaV.rs * veg[t] +
        # BetaW.rs * win[t] +
        EpsilonI.rsI[i] +
        EpsilonT.rsI[t]
    }
  }
  
  # #### Second attempt ####
  # for (i in 1:n){
  #   # likelihood
  #   rs[i] ~ dbern(rsI[id[i], year[i]])
  # 
  #   # constraints
  #   logit(rsI[id[i], year[i]]) <- logit(Mu.rsI[age[i]]) +
  #     EpsilonI.rsI[id[i]] +
  #     EpsilonT.rsI[year[i]]
  # }
  
  # use parameters estimated from individual data above
  # to predict age-specific reproductive success (rsA) here!
  for (a in 1:nAge){
    for (t in 1:nYear){
      logit(rsA[a, t]) <- logit(Mu.rsA[a]) + # rsA becomes s.PY in the population model!
        # BetaD.rs * dens[t] +
        # BetaV.rs * veg[t] +
        # BetaW.rs * win[t] +
        EpsilonT.rsA[t]
    }
  }
    
  ##### Priors ####
  # priors for fixed effects
  for (a in 1:nAge){
    Mu.rsI[a] ~ dunif(0, 1)
    Mu.rsA[a] ~ dunif(0, 1)
  }
  
  # Beta.dens ~ dunif(-2, 2) # could be dunif(-5, 5) if need be
  # Beta.veg  ~ dunif(-2, 2) # could be dunif(-5, 5) if need be
  # Beta.win  ~ dunif(-2, 2) # could be dunif(-5, 5) if need be
  
  # priors for random effects
  for (i in 1:nID.rs){
    EpsilonI.rsI[i] ~ dnorm(0, sd = SigmaI.rsI)
  }
  
  for (t in 1:nYear){
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
paraNimble <- function(seed, myCode, myConst, myData,
                       rs = rs, enData = enData, testRun){
  
  library(nimble)
  
  n = myConst$n
  nID.rs = myConst$nID.rs
  nYear = myConst$nYear
  nAge = myConst$nAge
  
  # assign initial values
  # TODO: UPDATE SIMULATE INITS FUNCTION!
  # myInits <- simulateInits(...)
  
  # assemble model
  myMod <- nimbleModel(code = myCode,
                       data = myData,
                       constants = myConst
                       # inits = myInits
                       )
  
  # select parameters to monitor
  params = c(# RS model
             "Mu.rsI", # "Mu.rsA"
             # "Beta.dens", "Beta.veg", "Beta.win",
             "EpsilonI.rsI", "EpsilonT.rsI",
             "SigmaI.rsI", "SigmaT.rsI"
  )
  
  # MCMC settings
  mySeed  <- 1:3
  nchains <- 3
  
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
  myMCMC <- buildMCMC(cModel, monitors = params, enableWAIC = T)
  cmyMCMC <- compileNimble(myMCMC, project = myMod)
  samples <- runMCMC(cmyMCMC,
                     niter = niter,
                     nburnin = nburnin,
                     thin = nthin,
                     setSeed = mySeed,
                     samplesAsCodaMCMC = T,
                     summary = T,
                     WAIC = T)
  
  return(samples)
}


## Run model -------------------------------------------------------------------

# serialized
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

