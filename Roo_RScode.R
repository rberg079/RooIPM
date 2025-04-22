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
source("wrangleData_env.R")
source("wrangleData_rs.R")

env <- wrangleData_env(dens.data = "data/WPNP_Methods_Results_January2025.xlsx",
                       veg.data  = "data/biomass data April 2009 - Jan 2025_updated Feb2025.xlsx",
                       wea.data  = "data/Prom_Weather_2008-2023_updated Jan2025 RB.xlsx",
                       wind.data = "data/POWER_Point_Daily_20080101_20241231_10M.csv")

rs <- wrangleData_rs(rs.data = "data/RSmainRB_Mar25.xlsx",
                     obs.data = "data/PromObs_2008-2019.xlsx",
                     prime = c(4:9), known.age = TRUE, cum.surv = TRUE, surv.sep1 = TRUE)

# check that all years are represented
setequal(1:17, unique(rs$year)) # should be TRUE!

# create Nimble lists
myData <-  list(y = rs$y,
                age = rs$age,
                dens = env$dens,
                veg = env$veg,
                win = env$win)

myConst <- list(N.id = rs$N.id,
                N.year = rs$N.year,
                N.age = rs$N.age)
# TODO: deal with missing environment


## Model -----------------------------------------------------------------------

myCode = nimbleCode({
  
  ##### Likelihood ####
  for (i in 1:N.id){
    for (t in 1:N.year){
      y[i, t] ~ dbern(p[i, t])
      logit(p[i, t]) <- mu.RS + B.age[age[i, t]] +
      B.dens * dens[t] + B.veg * veg[t] + B.win * win[t] +
      B.id[id[i]] + B.year[year[t]]
    }
  }
  
  for (i in 1:N.age){
    for (t in 1:N.year){
      y[a, t] ~ dbern(p[a, t])
      logit(p[a, t]) <- mu.RS + B.age[a] +
      B.dens * dens[t] + B.veg * veg[t] + B.win * win[t]
    }
  }
    
  ##### Priors ####
  # priors for fixed effects
  mu.RS  ~ dnorm(0, 1)
  B.age  ~ dnorm(0, 1)
  B.dens ~ dnorm(0, 1)
  B.veg  ~ dnorm(0, 1)
  B.win  ~ dnorm(0, 1)
  
  
  # priors for random effects
  for (j in 1:N.id){
    B.id[j] ~ dnorm(0, sigma[1])
  }
  
  for (t in 1:N.year){
    B.year[t] ~ dnorm(0, sigma[2])
  }
  
  # priors for sigma
  for (k in 1:2){
    sigma[i] ~ dunif(0, 100)
  }
  
})


## Assemble --------------------------------------------------------------------

# to serialize
# create Nimble function
paraNimble <- function(seed, myCode, myConst, myData,
                       rs = rs, env = env, testRun){
  
  library(nimble)
  
  N = myConst$N
  N.id = myConst$N.id
  N.year = myConst$N.year
  
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
             "mu.RS", "B.age", "B.dens", "B.veg", "B.win",
             "B.id", "B.year", "sigma"
  )
  
  # MCMC settings
  mySeed  <- 1
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
                     samplesAsCodaMCMC = T,
                     niter = niter,
                     nburnin = nburnin,
                     thin = nthin,
                     summary = T,
                     WAIC = T)
  
  return(samples)
}

