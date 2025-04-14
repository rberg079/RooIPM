# 8 April 2025
# Process model for roo IPM
# 1 census per year, on Sept 1 ish

## Set up ----------------------------------------------------------------------

# load packages
library(tidyverse)
library(lubridate)
library(beepr)
library(here)
# library(boot)
library(coda)
# library(foreach)
# library(doParallel)
# library(parallel)
library(nimble)
# registerDoParallel(3)

# load data
source("wrangleData_env.R")
source("wrangleData_surv.R")

env <- wrangleData_env(dens.data = "data/WPNP_Methods_Results_January2025.xlsx",
                       veg.data  = "data/biomass data April 2009 - Jan 2025_updated Feb2025.xlsx",
                       wea.data  = "data/Prom_Weather_2008-2023_updated Jan2025 RB.xlsx",
                       wind.data = "data/POWER_Point_Daily_20080101_20241231_10M.csv")

surv <- wrangleData_surv(surv.data = "data/PromSurvivalOct24.xlsx",
                         yafs.data = "data/RSmainRB_Mar25.xlsx")

# create Nimble lists
ntimes <- 17
nADs   <- 18
nAgeC  <- 5

myData  <- list(obs = surv$obs,
                state = surv$state,
                age = surv$age,
                ageC = surv$ageC,
                dens = env$dens,
                densE = env$densE,
                veg = env$veg,
                vegE = env$vegE)

myConst <- list(ntimes = ntimes,       # TODO: get from one of the wrangles?
                nind = surv$nind,
                nADs = nADs,           # TODO: get from wrangleData_rs
                nAgeC = surv$nAgeC,
                noAge = surv$noAge,
                nNoAge = surv$nNoAge,
                nNoVeg = env$nNoVeg,
                first = surv$first,
                last = surv$last,
                W = surv$W,
                DF = surv$DF)

# Switches/toggles
testRun <- TRUE # or FALSE


## Model -----------------------------------------------------------------------

myCode = nimbleCode({
  
  ## PROCESS MODEL
  ## ----------------------------------------
  
  for (t in 1:(ntimes-1)){
    YAF[t+1] ~ dbin(b[t] * s.PY[t], sum(AD[3:nADs,t])) 
    
    SA[1,t+1] ~ dbin(s.YAF[t], YAF[t])
    SA[2,t+1] ~ dbin(s.SA[1,t], SA[1,t])
    
    AD[1:2,t+1] <- 0
    AD[3,t+1] ~ dbin(s.SA[2,t], SA[2,t])
    
    for (a in 4:nADs){
      AD[a,t+1] ~ dbin(s.AD[a,t], AD[a,t])
    }
    Ntot[t+1] <- YAF[t+1] + sum(SA[1:2,t+1]) + sum(AD[3:nADs,t+1])
  }
  
  # priors
  for(t in 1:(ntimes-1)){
    b[t]        ~ dunif(0.5, 1)
    s.PY[t]     ~ dunif(0.1, 1)
    s.YAF[t]    <- s[1,t]
    s.SA[1,t]   <- s[2,t]
    s.SA[2,t]   <- s[2,t]
    
    for(a in 3:6){  # prime-aged
      s.AD[a,t] <- s[3,t]
    }
    
    for(a in 7:9){  # pre-senescent
      s.AD[a,t] <- s[4,t]
    }
    
    for(a in 10:nADs){  # senescent
      s.AD[a,t] <- s[5,t]
    }
  }
  
  ## CAPTURE-MARK-RECAPTURE MODEL (CJS)
  ## ----------------------------------------
  
  ## 1. Survival function
  for (a in 1:nAgeC){                               
    for (t in 1:(ntimes-1)){
      logit(s[a, t]) <- B.age[a] +
        B.dens[a] * dens.hat[t] +
        B.veg[a] * veg.hat[t] +
        # B.densVeg[a] * (dens.hat[t] * veg.hat[t]) +
        # B.vegRoo[a] * (veg.hat[t] / dens.hat[t]) +
        gamma[t, a]
    } # t
  } # a
  
  for (t in 1:ntimes){
    dens.hat[t] ~ dnorm(dens[t], sd = densE[t])
    veg.hat[t] ~ dnorm(veg[t], sd = vegE[t])
  }
  
  for (mt in 1:nNoVeg){
    veg[mt] ~ dnorm(0, 0.001)
  }
  
  # Estimate missing ages
  for (i in 1:nNoAge){
    ageM[i] ~ T(dnegbin(0.25,1.6), 3, 20)
    age[noAge[i], first[noAge[i]]] <- round(ageM[i]) + 1
    for (t in (first[noAge[i]]+1):ntimes){
      age[noAge[i], t] <- age[noAge[i], t-1] + 1
    } # t
  } # i
  
  # Priors for fixed effects
  for(a in 1:nAgeC){
    B.age[a] ~ dlogis(0, 1)
    B.dens[a] ~ dnorm(0, 0.001)
    B.veg[a] ~ dnorm(0, 0.001)
    # B.densVeg[a] ~ dnorm(0, 0.001)
  }
  
  # Variance-Covariance matrix
  for (i in 1:nAgeC){
    zero[i] <- 0
    xi[i] ~ dunif(0, 2)
  } # i
  
  for (t in 1:(ntimes-1)){
    eps.raw[t, 1:nAgeC]  ~ dmnorm(zero[1:nAgeC], Tau.raw[1:nAgeC,1:nAgeC])
    for (i in 1:nAgeC){
      gamma[t, i] <- xi[i] * eps.raw[t,i]
    } # i
  } # t
  
  # Priors for precision matrix
  Tau.raw[1:nAgeC, 1:nAgeC] ~ dwish(W[1:nAgeC, 1:nAgeC], DF)
  Sigma.raw[1:nAgeC, 1:nAgeC] <- inverse(Tau.raw[1:nAgeC, 1:nAgeC])
  
  # Uniform covariance matrix
  # for (a in 1:nAgeC){
  #   zero[a] <- 0
  #   sd.yr[a] ~ dunif(0, 5)
  #   cov.yr[a, a] <- sd.yr[a]*sd.yr[a]
  # }
  #
  # for(a in 1:(nAgeC-1)){
  #   for(a2 in (a+1):nAgeC){
  #     cor.yr[a, a2] ~ dunif(-1, 1)
  #     cov.yr[a2, a] <- sd.yr[a] * sd.yr[a2] * cor.yr[a, a2]
  #     cov.yr[a, a2] <- cov.yr[a2, a]
  #   }
  # } #i
  #
  # for (t in 1:(ntimes-1)) {
  #   gamma[t, 1:nAgeC]  ~ dmnorm(zero[1:nAgeC], cov = cov.yr[1:nAgeC, 1:nAgeC])
  # } #t
  
  
  ## 2. Observation
  for (t in 1:ntimes){
    for(i in 1:nind){
      logit(p[i, t]) <- mu.p + year.p[t]
    }
    year.p[t] ~ dnorm(0, sd = sd.p)
  }
  
  mu.p <- log(mean.p / (1-mean.p))
  mean.p ~ dunif(0, 1)
  sd.p ~ dunif(0, 10) # try dunif(0.01, 10) if sd.p causes trouble
  
  
  ### 3. Likelihood
  for (i in 1:nind){
    for (t in (first[i] + 1):last[i]){ # TODO: double-check that the first "first[i]" = first year in IPM
      # State process
      state[i, t] ~ dbern(mu1[i, t])
      mu1[i, t] <- s[ageC[age[i, t]], t-1] * state[i, t-1]
      
      # Observation process
      obs[i, t] ~ dbern(mu2[i, t])
      mu2[i, t] <- p[i, t] * state[i, t]
    } # t
  } # i
  
}) # nimbleCode


## Assemble --------------------------------------------------------------------

source("simulateInits.R")
myInits <- simulateInits(ntimes = ntimes, nADs = nADs, nAgeC = nAgeC,
                         dens = env$dens, veg = env$veg, nNoAge = surv$nNoAge)

# monitors
params = c(# CJS model
           'dens.hat', 'veg.hat', 'ageM',              # latent states
           'B.age', 'B.dens', 'B.veg', # 'B.densVeg',  # covariate effects
           'mean.p', 'year.p', 'sd.p',                 # observation parameters
           'gamma', 'xi', 'Sigma.raw',                 # correlated random effects
           # 'gamma', 'sd.yr', 'cor.yr'                # uniform random effects
           
           # Process model
           's', 'b', 's.PY', 's.YAF', 's.SA', 's.AD',  # yearly vita rates
           'YAF', 'SA', 'AD', 'Ntot')                  # population sizes


# # to run parallel chains
# # create Nimble function
# paraNimble <- function(seed, myCode, myConst, myData,
#                        n.burn = 1000, n.tin = 1){ # TODO: tweak MCMC params
#   
#   library(nimble)
#   
#   # assign initial values
#   source("simulateInits.R")
#   myInits <- simulateInits(ntimes = ntimes,
#                            nAgeC = nAgeC,
#                            nADs = nADs)
#   
#   # assemble model
#   myMod <- nimbleModel(code = myCode,
#                        data = myData,
#                        constants = myConst,
#                        inits = myInits())
#   
#   # select parameters to monitor
#   params = c(# CJS model
#              'dens.hat', 'veg.hat', 'ageM',              # latent states
#              'B.age', 'B.dens', 'B.veg', # 'B.densVeg',  # covariate effects
#              'mean.p', 'year.p', 'sd.p',                 # observation parameters
#              'gamma', 'xi', 'Sigma.raw',                 # correlated random effects
#              # 'gamma', 'sd.yr', 'cor.yr'                # uniform random effects
#              
#              # Process model
#              's', 'b', 's.PY', 's.YAF', 's.SA', 's.AD',  # yearly vita rates
#              'YAF', 'SA', 'AD', 'Ntot')                  # population sizes
#   
#   # select MCMC settings
#   cModel <- compileNimble(myMod)
#   myMCMC <- buildMCMC(cModel, monitors = vars, enableWAIC = T)
#   cmyMCMC <- compileNimble(myMCMC, project = myMod)
#   samples <- runMCMC(cmyMCMC,
#                      samplesAsCodaMCMC = T,
#                      niter = n.burn + 1000*n.tin,
#                      nburnin = n.burn,
#                      thin = n.tin,
#                      summary = T,
#                      WAIC = T)
#   
#   return(samples)
# }


## Run model -------------------------------------------------------------------

# MCMC specs
mySeed  <- 1
nthin   <- 1
nchains <- 1

if(testRun){
  niter   <- 10
  nburnin <- 0
}else{
  niter   <- 10000
  nburnin <- 2000
}

# run
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
                      setSeed = mySeed)


## Plots -----------------------------------------------------------------------

library(coda)
library(ggplot2)
library(MCMCvis)

# MCMC output
out.mcmc <- as.mcmc.list(samples)

summary(out.mcmc) # cannot handle NAs
MCMCsummary(out.mcmc, params = c('YAF', 'SA'), n.eff = TRUE, round = 2)

par(mfrow = c(4, 1))
plot(out.mcmc[, paste0('YAF[', 1:ntimes, ']')])
plot(out.mcmc[, paste0('SA[', rep(1:2, each = ntimes), ', ', rep(1:ntimes, times = 2), ']')])
# plot(out.mcmc[, paste0('AD[', rep(1:nADs, each = ntimes), ', ', rep(1:ntimes, times = nADs), ']')])

# assemble posterior samples
out.mat <- as.matrix(samples)

# parameters to include
table.params <- c(
  paste0('YAF[', 1:ntimes, ']'),
  paste0('SA[', rep(1:2, each = ntimes), ', ', rep(1:ntimes, times = 2), ']'),
  paste0('AD[', rep(1:nADs, each = ntimes), ', ', rep(1:ntimes, times = nADs), ']'))

# table.params <- list(
#   yaf = c(paste0('YAF[', 1:ntimes, ']')),
#   sa  = c(paste0('SA[', rep(1:2, each = ntimes), ', ', rep(1:ntimes, times = 2), ']')),
#   ad  = c(paste0('AD[', rep(1:nADs, each = ntimes), ', ', rep(1:ntimes, times = nADs), ']')))

# table with posterior summaries
post.table <- data.frame(Parameter = table.params, Estimate = NA)

for(i in 1:length(table.params)){
  est <- out.mat[, table.params[i]]
  post.table$Estimate[i] <- paste0(round(median(est, na.rm = T), digits = 2), ' [',
                                   round(quantile(est, 0.025, na.rm = T), digits = 2), ', ',
                                   round(quantile(est, 0.975, na.rm = T), digits = 2), ']')
}

# plot results
ntot <- grep("^Ntot\\[", colnames(samples)); ntot
yaf <- grep("^YAF\\[", colnames(samples)); yaf
sa <- grep("^SA\\[", colnames(samples)); sa
ad <- grep("^AD\\[", colnames(samples)); ad

var <- ntot

df <- data.frame(
  Year = 1:length(var),
  Mean = apply(out.mat[, var, drop = FALSE], 2, mean, na.rm = TRUE),
  Lower = apply(out.mat[, var, drop = FALSE], 2, quantile, probs = 0.025, na.rm = TRUE),
  Upper = apply(out.mat[, var, drop = FALSE], 2, quantile, probs = 0.975, na.rm = TRUE)
)

ggplot(df, aes(x = Year, y = Mean)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "skyblue", alpha = 0.4) +
  geom_line(color = "blue", linewidth = 1) +
  ylab("Parameter value") +
  xlab("Year") +
  theme_bw()

# ggsave("IPM.jpeg", scale = 1, width = 18.0, height = 9.0, units = c("cm"), dpi = 600)

