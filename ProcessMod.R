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
    s.YAF[t]    ~ dunif(0, 1)
    s.SA[1,t]   ~ dunif(0.5, 1)
    s.SA[2,t]   ~ dunif(0.5, 1)
    
    for(a in 3:6){ # prime-aged
      s.AD[a,t] ~ dunif(0.6, 1)
    }
    
    for(a in 7:9){ # pre-senescent
      s.AD[a,t] ~ dunif(0.5, 1)
    }
    
    for(a in 10:nADs){ # senescent
      s.AD[a,t] ~ dunif(0.1, 1)
    }
    
  }
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

