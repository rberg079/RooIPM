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
myData <-  list(rs    = rsData$survS1,
                id    = rsData$id,
                year  = rsData$year,
                age.R = rsData$age.R,
                dens  = enData$dens,
                veg   = enData$veg,
                win   = enData$win)

myConst <- list(nRS   = rsData$nRS,
                nID.R  = rsData$nID.R,
                nYear = rsData$nYear,
                nAge  = rsData$nAge,
                nAgeC = rsData$nAgeC)
# TODO: deal with missing environment

# Switches/toggles
testRun <- TRUE # or FALSE


## Model -----------------------------------------------------------------------

myCode = nimbleCode({
  
  #### Likelihood & constraints ####
  # individual rs function
  Mu.rsI[age[1]] <- 0
  Mu.rsI[age[2]] <- 0
  Mu.rsA[1] <- 0
  Mu.rsA[2] <- 0
  
  for(x in 1:nRS){
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
  for(a in 3:nAge){
    Mu.rsI[a] ~ dunif(0, 1)
    Mu.rsA[a] ~ dunif(0, 1)
  }
  
  # Beta.dens ~ dunif(-2, 2) # could be dunif(-5, 5) if need be
  # Beta.veg  ~ dunif(-2, 2) # could be dunif(-5, 5) if need be
  # Beta.win  ~ dunif(-2, 2) # could be dunif(-5, 5) if need be
  
  # priors for random effects
  for(i in 1:nID.R){
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

nchains   <- 3
seedMod   <- 1:nchains
seedInits <- 1

# assign initial values
source("simulateInits.R")
set.seed(seedInits)
myInits <- list()
for(c in 1:nchains){
  myInits[[c]] <- simulateInits(
    nRS = rsData$nRS,
    # nID.S = myConst$nID.S,
    nID.R = myConst$nID.R,
    nYear = myConst$nYear,
    nAge = myConst$nAge,
    nAgeC = myConst$nAgeC,
    
    age.R = myData$age.R,
    dens = myData$dens,
    veg = myData$veg,
    win = myData$win,
    
    # nNoAge = myConst$nNoAge,
    # nNoDens = myConst$nNoDens,
    nNoVeg = myConst$nNoVeg,
    nNoWin = myConst$nNoWin
  )
}

# select parameters to monitor
params = c("Mu.rsI", "Mu.rsA",
           # "Beta.dens", "Beta.veg", "Beta.win",
           "EpsilonI.rsI", "EpsilonT.rsI", "EpsilonT.rsA",
           "SigmaI.rsI", "SigmaT.rsI", "SigmaT.rsA"
)

# select MCMC settings
if(testRun){
  niter   <- 10
  nburnin <- 0
  nthin   <- 1
}else{
  niter   <- 10000
  nburnin <- 6000
  nthin   <- 1
}


## Run model -------------------------------------------------------------------

start <- Sys.time()
samples <- nimbleMCMC(code = myCode,
                      data = myData,
                      constants = myConst,
                      inits = myInits,
                      monitors = params,
                      niter = niter,
                      nburnin = nburnin,
                      nchains = nchains,
                      thin = nthin,
                      samplesAsCodaMCMC = T,
                      setSeed = seedMod)
dur <- Sys.time() - start; dur
beep(2)

library(MCMCvis)
out.mcmc <- samples %>% map(~as.mcmc(.$samples)) %>% as.mcmc.list()

MCMCsummary(out.mcmc,
            params = c('Mu.rsI', 'Mu.rsA', 'SigmaI.rsI', 'SigmaT.rsI', 'SigmaT.rsA'),
            n.eff = TRUE, round = 2)

