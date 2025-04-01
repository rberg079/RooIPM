# 31 March 2025
# Run survival model

## Set up ----------------------------------------------------------------------

# load packages
library(tidyverse)
library(lubridate)
library(beepr)
library(here)
library(boot)
library(coda)
library(foreach)
library(doParallel)
library(parallel)
library(nimble)
registerDoParallel(3)

# load data
source(here("PrepSURV.R"))
remove(surv, dens, densE, i, nNoAge, nNoDens, nNoVeg,
       noAge, noInfo, veg, vegE, fill_ages)

# # old data
# obs <- read_csv('data/obsF.csv') %>% as.matrix()
# state <- read_csv('data/stateF.csv') %>% as.matrix()
# age <- read_csv('data/ageF.csv') %>% as.matrix()
# id <- read_csv('data/idF.csv')

# limit to known-age inds for now!
uka <- id$uka

age   <- as.matrix(age[!uka,])
obs   <- as.matrix(obs[!uka,])
state <- as.matrix(state[!uka,])

first <- first[!uka]
last  <- last[!uka]
nind  <- nrow(state)

mydata  <- list(obs = obs, state = state, age = age, ageC = ageC)

myconst <- list(nind = nind, ntimes = ntimes, nAge = nAge, 
                first = first, last = last, W = diag(nAge), DF = nAge)


## Model -----------------------------------------------------------------------

myCode = nimbleCode({
  ##### 1. Survival ####
  # Survival function
  for (i in 1:nind){
    for (t in first[i]:(last[i]-1)){
      logit(s[i,t]) <- B.age[ageC[age[i,t]]] +
        gamma[t,ageC[age[i,t]]]
    } #t
  } #i
  
  # Priors for fixed effects
  for(a in 1:nAge){
    B.age[a] ~ dlogis(0,1)
  }
  
  # Variance-Covariance matrix
  for (i in 1:nAge){
    zero[i] <- 0
    xi[i] ~ dunif(0,2)
  } #i
  
  for (t in 1:(ntimes-1)){
    eps.raw[t,1:nAge]  ~ dmnorm(zero[1:nAge], Tau.raw[1:nAge,1:nAge])
    for (i in 1:nAge){
      gamma[t,i] <- xi[i] * eps.raw[t,i]
    } #i
  } #t
  
  # Priors for precision matrix
  Tau.raw[1:nAge, 1:nAge] ~ dwish(W[1:nAge,1:nAge], DF)
  Sigma.raw[1:nAge, 1:nAge] <- inverse(Tau.raw[1:nAge,1:nAge])
  
  
  ##### 2. Observation ####
  for (t in 1:ntimes){
    for(i in 1:nind){
      logit(p[i,t]) <- mu.p + year.p[t]
    }
    year.p[t] ~ dnorm(0, sd = sd.p)
  }
  
  mu.p <- log(mean.p / (1-mean.p))
  mean.p ~ dunif(0,1)
  sd.p ~ dunif(0,10)
  
  
  ##### 3. Likelihood ####
  for (i in 1:nind){
    for (t in (first[i] + 1):last[i]){
      # State process
      state[i,t] ~ dbern(mu1[i,t])
      mu1[i,t] <- s[i,t-1] * state[i,t-1]
      
      # Observation process
      obs[i,t] ~ dbern(mu2[i,t])
      mu2[i,t] <- p[i,t] * state[i,t]
    } #t
  } #i
  
})


## Assemble --------------------------------------------------------------------

# create Nimble function
paraNimble <- function(seed, myCode, myconst, mydata,
                       n.burn = 1000, n.tin = 1){
  
  library(nimble)
  library(coda)
  
  nind = myconst$ntimes
  ntimes = myconst$ntimes
  
  age = mydata$age
  nAge = myconst$nAge
  # noAge = myconst$noAge
  # nNoAge = myconst$nNoAge
  
  first = myconst$first
  last = myconst$last 
  
  # veg = mydata$veg
  # dens = mydata$dens
  # nNoVeg = myconst$nNoVeg
  
  # assign initial values
  myinits <- function(i){
    l = list(# ageM = sample(3:8, size = nNoAge, replace = T),
             
             B.age = rnorm(nAge,0,0.5), # rnorm(nAge,0,0.25)
             # B.veg = rnorm(nAge,0,0.5),
             # B.dens = rnorm(nAge,0,1),
             # B.densVeg = rnorm(nAge,0,1),
             
             # veg.hat = ifelse(is.na(veg),rnorm(length(veg),0,.1),veg),
             # dens.hat = ifelse(is.na(dens),rnorm(length(dens),0,.1),dens),
             
             mean.p = runif(1,0.6,1),
             year.p = rnorm(ntimes,0,0.2),
             sd.p = rnorm(1,0.2,0.1),
             
             xi = rnorm(nAge,1,0.1),
             eps.raw = matrix(rnorm((ntimes-1)*nAge,0,0.1),
                              ncol = nAge, nrow = (ntimes-1))
             
             # cor.yr = diag(nAge)+0.01,
             # sd.yr = runif(nAge,0,1)
    )
    Tau.raw = diag(nAge) + rnorm(nAge^2,0,0.1)
    l$Tau.raw = inverse((Tau.raw + t(Tau.raw))/2)
    return(l)
  }
  
  # assemble model
  myMod <- nimbleModel(code = myCode,
                       data = mydata,
                       constants = myconst,
                       inits = myinits())
  
  # select parameters to monitor
  vars = c('year.p', 'mean.p', 'sd.p', 'state', 
           'B.age', 'gamma', 'xi', 'Sigma.raw'
  )
  
  # select MCMC settings
  cModel <- compileNimble(myMod)
  mymcmc <- buildMCMC(cModel, monitors = vars, enableWAIC = T)
  CmyMCMC <- compileNimble(mymcmc, project = myMod)
  samples <- runMCMC(CmyMCMC,
                     samplesAsCodaMCMC = T,
                     niter = n.burn + 1000*n.tin,
                     nburnin = n.burn,
                     thin = n.tin,
                     summary = T,
                     WAIC = T)
  
  return(samples)
}


## Run model -------------------------------------------------------------------

# start.t <- Sys.time()
# this_cluster <- makeCluster(3)
# chain_output <- parLapply(cl = this_cluster,
#                           X = 1:3,
#                           fun = paraNimble,
#                           myCode = myCode,
#                           myconst = myconst,
#                           mydata = mydata,
#                           n.burn = 500,
#                           n.tin = 10)

start.t <- Sys.time()
chain_output <- paraNimble(seed = 123,
                           myCode = myCode,
                           myconst = myconst,
                           mydata = mydata,
                           n.burn = 10000,
                           n.tin = 10)

beep(sound = 2)
# stopCluster(this_cluster)
dur = now() - start.t
dur

# reformat output
codaSamp <- chain_output %>% map(~as.mcmc(.x$samples)) %>% as.mcmc.list()
codaSamp <- codaSamp[,!grepl('state',colnames(codaSamp[[1]]))]

# obtain WAIC value
waic <- chain_output %>% map_dbl(~.x$WAIC)
waic

# save output
fit1 <- list(model = myCode, codaSamp = codaSamp, waic = waic, dur = dur)
# write_rds(fit1, 'out_fit1.rds', compress = 'xz')

