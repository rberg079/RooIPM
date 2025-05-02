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
enData <- wrangleData_en(dens.data = "data/abundanceData_Proteus.csv",
                         veg.data  = "data/biomass data April 2009 - Jan 2025_updated Feb2025.xlsx",
                         wea.data  = "data/Prom_Weather_2008-2023_updated Jan2025 RB.xlsx",
                         wind.data = "data/POWER_Point_Daily_20080101_20241231_10M.csv",
                         obs.data  = "data/PromObs_2008-2019.xlsx",
                         list      = "data/PromlistAllOct24.xlsx")

source("wrangleData_rs.R")
rsData <- wrangleData_rs(rs.data = "data/RSmainRB_Mar25.xlsx",
                         obs.data = "data/PromObs_2008-2019.xlsx",
                         known.age = TRUE, cum.surv = FALSE)

# check that all years are represented
setequal(1:17, unique(rsData$year)) # should be TRUE

# create Nimble lists
myData <-  list(B      = rsData$B,
                R      = rsData$survS1,
                id.R   = rsData$id,
                year.R = rsData$year,
                age.R  = rsData$age.R,
                dens   = enData$dens,
                veg    = enData$veg,
                win    = enData$win)

myConst <- list(nR     = rsData$nR,
                nID.R  = rsData$nID.R,
                nYear  = rsData$nYear,
                nAge   = rsData$nAge,
                nAgeC  = rsData$nAgeC)
# TODO: deal with missing environment

# Switches/toggles
testRun <- TRUE # or FALSE


## Model -----------------------------------------------------------------------

myCode = nimbleCode({
  
  #### Likelihood & constraints ####
  # yearly birth rate
  for(t in 1:(nYear-1)){
    B[t] ~ dbern(Bt[t])
    logit(Bt[t]) <- logit(Mu.Bt) + EpsilonT.Bt[t]
  }
  
  # individual RS function
  Mu.Ri[1] <- 0
  for(x in 1:nR){
    R[x] ~ dbern(Ri[x])
    logit(Ri[x]) <- logit(Mu.Ri[age.R[x]]) +
      # BetaD.R * dens[year.R[x]] +
      # BetaV.R * veg[year.R[x]] +
      # BetaW.R * win[year.R[x]] +
      EpsilonI.Ri[id.R[x]] +
      EpsilonT.Ri[year.R[x]]
  }
  
  # age-specific RS function
  # use parameters estimated from individual data above
  # to predict age-specific reproductive success (Ra) here!
  Mu.Ra[1] <- 0
  for(a in 1:nAge){
    for(t in 1:(nYear-1)){
      logit(Ra[a, t]) <- logit(Mu.Ra[a]) + # Ra used in Pop model
        # BetaD.R * dens[t] +
        # BetaV.R * veg[t] +
        # BetaW.R * win[t] +
        EpsilonT.Ra[t]
    }
  }
  
  ##### Priors ####
  # priors for fixed effects
  for(a in 2:nAge){
    Mu.Ri[a] ~ dunif(0, 1)
    Mu.Ra[a] ~ dunif(0, 1)
  }
  
  # BetaD.R ~ dunif(-2, 2) # could be dunif(-5, 5) if need be
  # BetaV.R  ~ dunif(-2, 2) # could be dunif(-5, 5) if need be
  # BetaW.R  ~ dunif(-2, 2) # could be dunif(-5, 5) if need be
  
  # priors for random effects
  for(i in 1:nID.R){
    EpsilonI.Ri[i] ~ dnorm(0, sd = SigmaI.Ri)
  }
  
  for(t in 1:(nYear-1)){
    EpsilonT.Ri[t] ~ dnorm(0, sd = SigmaT.Ri)
    EpsilonT.Ra[t] ~ dnorm(0, sd = SigmaT.Ra)
    EpsilonT.Bt[t] ~ dnorm(0, sd = SigmaT.Bt)
  }
  
  # priors for sigma
  SigmaI.Ri ~ dunif(0, 100)
  SigmaT.Ri ~ dunif(0, 100)
  SigmaT.Ra ~ dunif(0, 100)
  SigmaT.Bt ~ dunif(0, 100)
  
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
    nR = rsData$nR,
    nID.R = myConst$nID.R,
    nYear = myConst$nYear,
    nAge = myConst$nAge,
    nAgeC = myConst$nAgeC,
    id.R = myData$id.R,
    year.R = myData$year.R,
    age.R = myData$age.R
    # dens = myData$dens,
    # veg = myData$veg,
    # win = myData$win,
    # nNoDens = myConst$nNoDens,
    # nNoVeg = myConst$nNoVeg,
    # nNoWin = myConst$nNoWin
  )
}

# select parameters to monitor
params = c("Mu.Bt", "Mu.Ri", "Mu.Ra",
           # "BetaD.R", "BetaV.R", "BetaW.R",
           "EpsilonI.Ri", "EpsilonT.Ri", "EpsilonT.Ra", "EpsilonT.Bt",
           "SigmaI.Ri", "SigmaT.Ri", "SigmaT.Ra", "SigmaT.Bt"
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
out.mcmc <- as.mcmc.list(samples)

MCMCsummary(out.mcmc,
            params = c('Mu.Ri', 'Mu.Ra', 'SigmaI.Ri', 'SigmaT.Ri', 'SigmaT.Ra'),
            n.eff = TRUE, round = 2)

