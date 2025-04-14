# 8 April 2025
# Process model for roo IPM
# 1 census per year, on Sept 1 ish

## Set up ----------------------------------------------------------------------

# load packages
library(tidyverse)
library(nimble)

# Switches/toggles
testRun <- TRUE
#testRun <- FALSE

# create Nimble lists
ntimes <- 20
nADs   <- 18

mydata  <- list()
myconst <- list(ntimes = ntimes, nADs = nADs)


## Model -----------------------------------------------------------------------

myCode = nimbleCode({
  # process model
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
    
    for(a in 3:6){ # prime-aged
      s.AD[a,t] <- s[3,t]
    }
    
    for(a in 7:9){ # pre-senescent
      s.AD[a,t] <- s[4,t]
    }
    
    for(a in 10:nADs){ # senescent
      s.AD[a,t] <- s[5,t]
    }
    
  }
  
  #------------------------------------#
  # CAPTURE-MARK RECAPTURE MODEL (CJS) #
  #------------------------------------#
  
  ##### 1. Survival ####
  # Survival function
  for (a in 1:nAge_CJS){                               
    for (t in 1:(ntimes-1)){
      logit(s[a,t]) <- B.age[a] +
        dens.hat[t]*B.dens[a] +
        veg.hat[t]*B.veg[a] +
        # (dens.hat[t]*veg.hat[t])*B.densVeg[a] +
        # (veg.hat[t]/dens.hat[t])*B.vegRoo[a] +
        gamma[t,a]
    } #t
  } #a
  
  for (t in 1:ntimes){
    dens.hat[t] ~ dnorm(dens[t], sd = densE[t])
    veg.hat[t] ~ dnorm(veg[t], sd = vegE[t])
  }
  
  for (mt in 1:nNoVeg){
    veg[mt] <- dnorm(0,0.001)
  }
  
  # Estimate missing ages
  for (i in 1:nNoAge){
    ageM[i] ~ T(dnegbin(0.25,1.6),3,20)
    age[noAge[i],first[noAge[i]]] <- round(ageM[i])+1
    for (t in (first[noAge[i]]+1):ntimes){
      age[noAge[i],t] <- age[noAge[i],t-1]+1
    } #t
  } #i
  
  # Priors for fixed effects
  for(a in 1:nAge){
    B.age[a] ~ dlogis(0,1)
    B.dens[a] ~ dnorm(0,0.001)
    B.veg[a] ~ dnorm(0,0.001)
    # B.densVeg[a] ~ dnorm(0,0.001)
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
  
  # Uniform covariance matrix
  # for (a in 1:nAge){
  #   zero[a] <- 0
  #   sd.yr[a] ~ dunif(0,5)
  #   cov.yr[a,a] <- sd.yr[a]*sd.yr[a]
  # }
  #
  # for(a in 1:(nAge-1)){
  #   for(a2 in (a+1):nAge){
  #     cor.yr[a,a2] ~ dunif(-1,1)
  #     cov.yr[a2,a] <- sd.yr[a] * sd.yr[a2] * cor.yr[a,a2]
  #     cov.yr[a,a2] <- cov.yr[a2,a]
  #   }
  # } #i
  #
  # for (t in 1:(ntimes-1)) {
  #   gamma[t,1:nAge]  ~ dmnorm(zero[1:nAge], cov=cov.yr[1:nAge, 1:nAge])
  # } #t
  
  
  ##### 2. Observation ####
  for (t in 1:ntimes){
    for(i in 1:nind){
      logit(p[i,t]) <- mu.p + year.p[t]
    }
    year.p[t] ~ dnorm(0, sd = sd.p)
  }
  
  mu.p <- log(mean.p / (1-mean.p))
  mean.p ~ dunif(0,1)
  sd.p ~ dunif(0,10) # try dunif(0.01,10) if sd.p causes trouble
  
  
  ##### 3. Likelihood ####
  for (i in 1:nind){
    for (t in (first[i] + 1):last[i]){ #TODO: Double-check that the first "first[i]" = first year in IPM
      # State process
      state[i,t] ~ dbern(mu1[i,t])
      mu1[i,t] <- s[ageC[age[i,t]],t-1] * state[i,t-1]
      
      # Observation process
      obs[i,t] ~ dbern(mu2[i,t])
      mu2[i,t] <- p[i,t] * state[i,t]
    } #t
  } #i
  
}) # nimbleCode


## Assemble --------------------------------------------------------------------

source("simulateInits.R")
myinits <- simulateInits(ntimes = ntimes, nADs = nADs)

# monitors
params = c('YAF', 'SA', 'AD', 'Ntot',
           'b', 's.PY', 's.YAF', 's.SA', 's.AD')


## Run model -------------------------------------------------------------------

# MCMC specs
mySeed  <- 1

if(!testRun){
  niter   <- 2000
  nburnin <- 1000
}else{
  niter   <- 10
  nburnin <- 0
}

nthin   <- 1
nchains <- 1

# run
samples <- nimbleMCMC(code = myCode,
                      data = mydata,
                      constants = myconst,
                      inits = myinits,
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

par(mfrow = c(4,1))
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

