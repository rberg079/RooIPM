# 9 April 2025
# Run survival model

## Set up ----------------------------------------------------------------------

# load packages
library(tidyverse)
library(lubridate)
library(beepr)
library(here)
library(boot)
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

source("wrangleData_sv.R")
svData <- wrangleData_sv(surv.data = "data/PromSurvivalOct24.xlsx",
                         yafs.data = "data/RSmainRB_Mar25.xlsx")

# create Nimble lists
myData  <- list(obs = svData$obs,
                state = svData$state,
                age = svData$age,
                ageC = svData$ageC,
                dens = enData$dens,
                densE = enData$densE,
                veg = enData$veg,
                vegE = enData$vegE)

myConst <- list(nID.S = svData$nID.S,
                nYear = svData$nYear,
                nAgeC = svData$nAgeC,
                noAge = svData$noAge,
                nNoAge = svData$nNoAge,
                nNoVeg = enData$nNoVeg,
                first = svData$first,
                last = svData$last,
                W = svData$W,
                DF = svData$DF)

# # checks
# sapply(1:nrow(myData$state), function (i) myData$state[i, myConst$first[i]]) %>% table(useNA = 'a') # should all be 1
# sapply(1:nrow(myData$state), function (i) myData$state[i, myConst$last[i]]) %>% table(useNA = 'a')  # should be mostly 0s & NAs
# sapply(1:nrow(myData$state), function (i) myData$age[i, myConst$first[i]]) %>% table(useNA = 'a')   # should all be >= 1
# table(myConst$first >= myConst$last, useNA = 'a')                                                   # should all be F

# Switches/toggles
testRun <- TRUE # or FALSE


## Model -----------------------------------------------------------------------

myCode = nimbleCode({
  
  #### Likelihood ####
  for(i in 1:nID.S){
    for(t in (first[i] + 1):last[i]){
      # state process
      state[i, t] ~ dbern(Mu.Sp[i, t])
      Mu.Sp[i, t] <- sv[i, t-1] * state[i, t-1]
      # observation process
      obs[i, t] ~ dbern(Mu.Op[i, t])
      Mu.Op[i, t] <- O[i, t] * state[i, t]
    }
  }
  
  #### Constraints ####
  # survival function
  for(i in 1:nID.S){                               
    for(t in first[i]:(last[i]-1)){
      logit(sv[i, t]) <- BetaA.sv[ageC[age[i, t]]] +
        BetaD.sv[ageC[age[i, t]]] * dens.hat[t] +
        BetaV.sv[ageC[age[i, t]]] * veg.hat[t] +
        # BetaDV.sv[ageC[age[i, t]]] * (dens.hat[t] * veg.hat[t]) +
        # BetaVR.sv[ageC[age[i, t]]] * (veg.hat[t] / dens.hat[t]) +
        Gamma.sv[t, ageC[age[i,t]]]
    }
  }
  
  for(t in 1:nYear){
    dens.hat[t] ~ dnorm(dens[t], sd = densE[t])
    veg.hat[t] ~ dnorm(veg[t], sd = vegE[t])
  }
  
  for(m in 1:nNoVeg){
    veg[m] <- 0
    # veg[m] ~ dnorm(0, 0.001)
  }
  
  # estimate missing ages
  for(i in 1:nNoAge){
    ageM[i] ~ T(dnegbin(0.25, 1.6), 3, 20)
    age[noAge[i], first[noAge[i]]] <- round(ageM[i]) + 1
    for(t in (first[noAge[i]]+1):nYear){
      age[noAge[i], t] <- age[noAge[i], t-1] + 1
    }
  }
  
  # observation function
  for(t in 1:nYear){
    for(i in 1:nID.S){
      logit(O[i, t]) <- logit(Mu.O) + Epsilon.O[t]
    }
    Epsilon.O[t] ~ dnorm(0, sd = Sigma.O)
  }
  
  #### Priors ####
  # for fixed effects
  for(a in 1:nAgeC){
    BetaA.sv[a] ~ dnorm(0, sd = 2) # or dlogis(0, 1)
    BetaD.sv[a] ~ dnorm(0, 0.001)
    BetaV.sv[a] ~ dnorm(0, 0.001)
    # BetaDV.sv[a] ~ dnorm(0, 0.001)
  }
  
  # for random effects
  # variance-covariance matrix
  for(i in 1:nAgeC){
    zero[i] <- 0
    Xi.sv[i] ~ dunif(0, 2)
  }
  
  for(t in 1:(nYear-1)){
    Epsilon.sv[t, 1:nAgeC] ~ dmnorm(zero[1:nAgeC], Tau.sv[1:nAgeC, 1:nAgeC])
    for(i in 1:nAgeC){
      Gamma.sv[t, i] <- Xi.sv[i] * Epsilon.sv[t,i]
    }
  }
  
  # precision matrix
  Tau.sv[1:nAgeC, 1:nAgeC] ~ dwish(W[1:nAgeC, 1:nAgeC], DF)
  Sigma.sv[1:nAgeC, 1:nAgeC] <- inverse(Tau.sv[1:nAgeC, 1:nAgeC])
  
  # observation
  Mu.O ~ dunif(0.01, 0.99) # or dunif(0, 1)
  Sigma.O ~ dunif(0.01, 10) # or dunif(0, 10)
  
})


## Assemble --------------------------------------------------------------------

nchains   <- 3
seedMod   <- 1:nchains
seedInits <- 1

# assign initial values
# source("simulateInits.R")
# set.seed(seedInits)
myInits <- list()
for(c in 1:nchains){
  myInits[[c]] <- list(
    ageM = sample(3:8, size = myConst$nNoAge, replace = T),
    
    BetaA.sv = rnorm(myConst$nAgeC, 0, 1),
    BetaV.sv = rnorm(myConst$nAgeC, 0, 1),
    BetaD.sv = rnorm(myConst$nAgeC, 0, 1),
    # BetaDV.sv = rnorm(myConst$nAgeC, 0, 1),
    
    veg.hat = ifelse(is.na(myData$veg), rnorm(length(myData$veg), 0, .1), myData$veg),
    dens.hat = ifelse(is.na(myData$dens), rnorm(length(myData$dens), 0, .1), myData$dens),
    
    O <- matrix(runif(myConst$nID.S * myConst$nYear, 0.1, 0.9),
                 nrow = myConst$nID.S, ncol = myConst$nYear),
    
    Mu.O = runif(1, 0.1, 0.9),
    Epsilon.O = rnorm(myConst$nYear, 0, 0.2),
    Sigma.O = runif(1, 0.01, 2), # or rnorm(1, 0.2, 0.1)
    
    Xi.sv = rnorm(myConst$nAgeC, 1, 0.1),
    Epsilon.sv = matrix(rnorm((myConst$nYear-1)*myConst$nAgeC, 0, 0.1),
                        ncol = myConst$nAgeC, nrow = myConst$nYear-1),
    
    Tau.sv = diag(myConst$nAgeC) + rnorm(myConst$nAgeC^2, 0, 0.1),
    Tau.sv = inverse((Tau.sv + t(Tau.sv))/2) # may be a typo?
  )
}

# select parameters to monitor
params = c('Epsilon.O', 'Sigma.O', 'state', 'ageM',
           'BetaA.sv', 'BetaD.sv', 'BetaV.sv', # 'BetaDV.sv',
           'dens.hat', 'veg.hat', 'densE', 'vegE',
           'Gamma.sv', 'Xi.sv', 'Sigma.sv')

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

# reformat output
out.mcmc <- samples %>% map(~ as.mcmc(.x$samples)) %>% as.mcmc.list()
out.mcmc <- out.mcmc[, !grepl('state', colnames(out.mcmc[[1]]))]

# obtain WAIC value
waic <- samples %>% map_dbl(~.x$WAIC$WAIC)  # serialized
waic

# save output
fit <- list(model = myCode, out.mcmc = out.mcmc, waic = waic, dur = dur)
# write_rds(fit, 'results/out_fit.rds', compress = 'xz')


## Checks ----------------------------------------------------------------------

summary(out.mcmc) # cannot handle any NAs

library(MCMCvis)
library(corrplot)
MCMCsummary(out.mcmc, params = c('BetaA.sv', 'BetaD.sv', 'BetaV.sv'), n.eff = TRUE, round = 3)
MCMCsummary(out.mcmc, params = c('Epsilon.O','Sigma.O'), n.eff = TRUE, round = 3)
MCMCsummary(out.mcmc, params = c('Sigma.sv'), n.eff = TRUE, round = 3)
MCMCsummary(out.mcmc, params = c('ageM'), n.eff = TRUE, round = 3)

# assess MCMC convergence
# ...visually
nYear <- myConst$nYear
nAgeC <- myConst$nAgeC

par(mar = c(1,1,1,1))
plot(out.mcmc[, paste0('BetaA.sv[',1:nAgeC,']')])
plot(out.mcmc[, paste0('BetaV.sv[',1:nAgeC,']')])
plot(out.mcmc[, paste0('BetaD.sv[',1:nAgeC,']')])
plot(out.mcmc[, paste0('BetaDV.sv[',1:nAgeC,']')])

plot(out.mcmc[, paste0('Epsilon.O[',1:nYear,']')])
plot(out.mcmc[, 'Sigma.O'])

plot(out.mcmc[, paste0('Sigma.sv[',1:nAgeC,', ',1:nAgeC,']')])
plot(out.mcmc[, paste0('ageM[',1:nNoAge,']')])

# # ...formally
# gelman.diag(out.mcmc[, paste0('BetaA.sv[',1:nAgeC,']')])
# gelman.diag(out.mcmc[, paste0('BetaV.sv[',1:nAgeC,']')])
# gelman.diag(out.mcmc[, paste0('BetaD.sv[',1:nAgeC,']')])
# gelman.diag(out.mcmc[, paste0('BetaDV.sv[',1:nAgeC,']')])
# 
# gelman.diag(out.mcmc[, paste0('Epsilon.O[',1:nYear,']')])
# gelman.diag(out.mcmc[, 'Sigma.O'])
# 
# gelman.diag(out.mcmc[, paste0('Sigma.sv[',1:nAgeC,', ',1:nAgeC,']')])
# gelman.diag(out.mcmc[, paste0('ageM[',1:nNoAge,']')])
# 
# # check Neff
# effectiveSize(out.mcmc[, paste0('BetaA.sv[',1:nAgeC,']')])
# effectiveSize(out.mcmc[, paste0('BetaV.sv[',1:nAgeC,']')])
# effectiveSize(out.mcmc[, paste0('BetaD.sv[',1:nAgeC,']')])
# effectiveSize(out.mcmc[, paste0('BetaDV.sv[',1:nAgeC,']')])
# 
# effectiveSize(out.mcmc[, paste0('Epsilon.O[',1:nYear,']')])
# effectiveSize(out.mcmc[, 'Sigma.O'])
# 
# effectiveSize(out.mcmc[, paste0('Sigma.sv[',1:nAgeC,', ',1:nAgeC,']')])
# effectiveSize(out.mcmc[, paste0('ageM[',1:nNoAge,']')])


## Plots -----------------------------------------------------------------------

betaz = out.mcmc %>% map(as.data.frame) %>% bind_rows()

# check for correlations among fixed effects
par(mfrow = c(1,1))
corrplot(cor(out.mcmc[, grepl('B.', colnames(out.mcmc[[1]]))] %>%
               map(as.data.frame) %>%
               bind_rows()
             , use = 'p'))

# check random effects among demographic rates
# check variance-correlation matrix, with Sigma.sv on diagonal
varCorrMatrix <- array(NA, dim = c(myConst$nAgeC, myConst$nAgeC, nrow(betaz)))

for(i in 1:myConst$nAgeC){
  varCorrMatrix[i,i,] <- betaz[, paste0('Xi.sv[', i,']')]*
    sqrt(betaz[, paste0('Sigma.sv[', i,', ', i,']')])
}

for(j in 1:(myConst$nAgeC-1)){
  for(i in (j+1):myConst$nAgeC){
    varCorrMatrix[j,i,] <- (betaz[, paste0('Sigma.sv[',i,', ',j,']')])/
      sqrt(betaz[, paste0('Sigma.sv[', j,', ', j,']')]*
             betaz[, paste0('Sigma.sv[', i,', ', i,']')])
  }
}

round(apply(varCorrMatrix, 1:2, mean, na.rm = T), 2)
round(apply(varCorrMatrix, 1:2, quantile, prob = 0.025, na.rm = T), 2)
round(apply(varCorrMatrix, 1:2, quantile, prob = 0.975, na.rm = T), 2)

# calculate survival probabilities
df = expand.grid(age = 1:22, year = 1:16)
sv.pred = matrix(NA, nrow = nrow(df), ncol = nrow(betaz))

for(i in 1:nrow(df)){
  sv.pred[i,] <- betaz[, paste0('BetaA.sv[', myData$ageC[df$age[i]], ']')] +
    betaz[,paste0('Gamma.sv[', df$year[i], ', ', myData$ageC[df$age[i]], ']')]
  df$ageC[i] <- myData$ageC[df$age[i]]
}

df$sv = inv.logit(apply(sv.pred, 1, mean))
df$sv.cil = inv.logit(apply(sv.pred, 1, quantile, 0.025))
df$sv.cih = inv.logit(apply(sv.pred, 1, quantile, 0.975))

# plot main results
df <- df %>% mutate(ageC = as.factor(ageC))

library(ggplot2)
library(scales)
plot <- ggplot(df, aes(x = year, y = sv)) +
  geom_ribbon(aes(fill = ageC, ymin = sv.cil, ymax = sv.cih), alpha = 0.2) +
  geom_path(aes(colour = ageC), show.legend = F) +
  labs(x = "Year", y = "Survival", fill = "Age class") +
  # scale_x_continuous(breaks = pretty_breaks()) +
  scale_y_continuous(breaks = pretty_breaks(),
                     limits = c(0,1)) +
  theme_bw(); plot

