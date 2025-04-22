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
year <- as.numeric(factor(rs$year)); year
setequal(1:17, unique(year)) # should be TRUE!

# create Nimble lists
myData <-  list(y = rs$survS1,
                id = factor(rs$id),
                year = year,
                age = rs$age,
                # prs = factor(rs$prs),
                dens = env$dens,
                veg = env$veg,
                win = env$win)

myConst <- list(N = length(id),
                N.id = max(id),
                N.year = max(year))


## Model -----------------------------------------------------------------------

myCode = nimbleCode({
  
  ##### 1. Reproductive success ####
  
  for (i in 1:N){
    ##### Likelihood ####
    # RS function
    y[i] ~ dbern(p[i])
    logit(p[i]) <- int + b.age*age[i] + # b.prs[prs[i]]
      b.veg*veg[i] + b.dens*dens[i] + b.win*win[i] +
      b.yr[year[i]] + b.id[id[i]]
    
    # missing values
    dens[i] ~ dnorm(0, 1)
    veg[i] ~ dnorm(0, 1)
    win[i] ~ dnorm(0, 1)
    
    # prs[i] ~ dcat(p.prs[1:4])
    
    ##### Priors ####
    # priors for fixed effects
    b.age ~ dnorm(0, 1)

    # b.prs[1] <- 0
    # b.prs[2] ~ dnorm(0, 1)
    # b.prs[3] ~ dnorm(0, 1)
    # b.prs[4] ~ dnorm(0, 1)    
    # p.prs[1:4] ~ ddirch(rep(0.25, 4))
    
    b.dens ~ dnorm(0, 1)
    b.veg ~ dnorm(0, 1)
    b.win ~ dnorm(0, 1)
    
    int ~ dnorm(0, 1)
    
    # priors for random effects
    for (j in 1:N.id){
      b.id[j] ~ dnorm(0, sigma[1])
    }
    
    for (t in 1:N.year){
      b.yr[t] ~ dnorm(0, sigma[2])
    }
    
    # priors for sigma
    for (k in 1:2){
      sigma[i] ~ dunif(0, 100)
    }
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
    # TODO!
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

















# write model to file
write(mod, "modFULLsurv21_cond.txt")

# fit model to DATA
inits <- function() {
  list(cuts = sort(rnorm(5, 0, 1)))
}

mod <- jags(model = "modFULLsurv21_cond.txt",
            data = list(y = y, id = id, year = year,
                        N = N, N.id = N.id, N.year = N.year, age = age,
                        leg = leg, teeth = teeth, prs = prs, cond = cond, # mass = mass,
                        xmed = xmed, veg = veg, dens = dens, win = win, # bday = bday,
                        veg2 = veg2, dens2 = dens2, win2 = win2,
                        year.mon = year.mon, firstAge = firstAge
                        ),
            param = c(
              "r", "cuts", "int.tt", "b.age.tt",
              "mu.cond", "int.cd", "b.age.cd", "b.prs.cd", # "b.leg.cd", 
              "b.xmed.cd", "b.veg.cd", "b.dens.cd", "b.win.cd", # "b.bday.cd",
              "p", "int.rs", "b.cond.rs", "b.age.rs", "b.prs.rs", "b.leg.rs", 
              "b.xmed.rs", "b.veg.rs", "b.dens.rs", "b.win.rs", # "b.bday.rs",
              "sigma", "pred.cond", "pred.rs"
              ),
            n.chains = 3,
            n.iter = 60000,
            n.burnin = 40000,
            parallel = T)
beep(sound = 2)

# save
save(file = "Results/modFULLsurv21_cond.Rdata", list = "mod")

# # load
# load("Results/modFULLsurv21_cond.Rdata")

# summarise model output
MCMCsummary(mod,
            params = c(
              "r", "cuts", "int.tt", "b.age.tt",
              "int.cd", "b.age.cd", "b.prs.cd", # "mu.cond", "b.leg.cd", 
              "b.xmed.cd", "b.veg.cd", "b.dens.cd", "b.win.cd", # "b.bday.cd",
              "int.rs", "b.cond.rs", "b.age.rs", "b.leg.rs", "b.prs.rs", # "p", 
              "b.xmed.rs", "b.veg.rs", "b.dens.rs", "b.win.rs", # "b.bday.rs",
              "sigma" # "pred.cond", "pred.rs"
              ),
            round = 2)

MCMCtrace(mod, 
          params = c(
            "r", "cuts", "int.tt", "b.age.tt",
            "int.cd", "b.age.cd", "b.prs.cd", # "mu.cond", "b.leg.cd", 
            "b.xmed.cd", "b.veg.cd", "b.dens.cd", "b.win.cd", # "b.bday.cd",
            "int.rs", "b.cond.rs", "b.age.rs", "b.leg.rs", "b.prs.rs", # "p", 
            "b.xmed.rs", "b.veg.rs", "b.dens.rs", "b.win.rs", # "b.bday.rs",
            "sigma" # "pred.cond", "pred.rs"
          ),
          pdf = FALSE)

par(mfrow = c(1,1))
MCMCplot(mod, 
         params = c(
           "r", "cuts", "int.tt", "b.age.tt",
           "int.cd", "b.age.cd", "b.prs.cd", # "mu.cond", "b.leg.cd", 
           "b.xmed.cd", "b.veg.cd", "b.dens.cd", "b.win.cd", # "b.bday.cd",
           "int.rs", "b.cond.rs", "b.age.rs", "b.leg.rs", "b.prs.rs", # "p", 
           "b.xmed.rs", "b.veg.rs", "b.dens.rs", "b.win.rs", # "b.bday.rs",
           "sigma" # "pred.cond", "pred.rs"
         ),
         ci = c(50, 90))

whiskerplot(mod,
            parameters = c("b.age.cd", "b.prs.cd", "b.xmed.cd",
                           "b.veg.cd", "b.dens.cd", "b.win.cd")) # "b.bday.ma"

whiskerplot(mod,
            parameters = c("b.cond.rs", "b.age.rs", "b.leg.rs", "b.prs.rs", "b.xmed.rs",
                           "b.veg.rs", "b.dens.rs", "b.win.rs")) # "b.bday.rs"

