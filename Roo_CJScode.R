# 9 April 2025
# Run survival model

## Set up ----------------------------------------------------------------------

# set toggles
testRun <- FALSE
parallelRun <- TRUE

# load packages
library(tidyverse)
library(lubridate)
library(beepr)
library(here)
library(boot)
library(coda)
library(nimble)
library(parallel)

# load data
source('wrangleData_en.R')
enData <- wrangleData_en(dens.data = "data/abundanceData_Proteus.csv",
                         veg.data  = "data/biomass data April 2009 - Jan 2025_updated Feb2025.xlsx",
                         wea.data  = "data/Prom_Weather_2008-2023_updated Jan2025 RB.xlsx",
                         wind.data = "data/POWER_Point_Daily_20080101_20241231_10M.csv",
                         obs.data  = "data/PromObs_2008-2019.xlsx",
                         list      = "data/PromlistAllOct24.xlsx")

source('wrangleData_sv.R')
svData <- wrangleData_sv(surv.data = "data/PromSurvivalOct24.xlsx",
                         yafs.data = "data/RSmainRB_Mar25.xlsx",
                         known.age = TRUE)

# to play around with age classes!
# ageC <- c(1,2,2,3,3,3,3,4,4,4, rep(5,30)) # default
ageC <- c(seq(from = 1, to = 20, by = 1), rep(20, times = 20)); ageC

# create Nimble lists
myData  <- list(obs   = svData$obs,
                state = svData$state,
                age.S = svData$age.S,
                ageC  = ageC,
                
                dens  = enData$dens,
                densE = enData$densE,
                veg   = enData$veg,
                vegE  = enData$vegE)

myConst <- list(nID.S = svData$nID.S,
                nYear = svData$nYear,
                nAgeC = max(ageC),
                first = svData$first,
                last  = svData$last,
                W     = diag(max(ageC)),
                DF    = max(ageC))

# # checks
# sapply(1:nrow(myData$state), function (i) myData$state[i, myConst$first[i]]) %>% table(useNA = 'a') # should all be 1
# sapply(1:nrow(myData$state), function (i) myData$state[i, myConst$last[i]]) %>% table(useNA = 'a')  # should be mostly 0s & NAs
# sapply(1:nrow(myData$state), function (i) myData$age[i, myConst$first[i]]) %>% table(useNA = 'a')   # should all be >= 1
# table(myConst$first >= myConst$last, useNA = 'a')                                                   # should all be F


## Model -----------------------------------------------------------------------

myCode = nimbleCode({
  
  #### Likelihood ####
  for(i in 1:nID.S){
    for(t in (first[i] + 1):last[i]){
      # state process
      state[i, t] ~ dbern(Mu.Sp[i, t])
      Mu.Sp[i, t] <- S[ageC[age.S[i, t-1]], t-1] * state[i, t-1]
      
      # observation process
      obs[i, t] ~ dbern(Mu.Op[i, t])
      Mu.Op[i, t] <- O[t] * state[i, t]
    }
  }
  
  #### Constraints ####
  # survival function
  for(a in 1:nAgeC){                               
    for(t in 1:(nYear-1)){
      logit(S[a, t]) <- BetaA.S[a] +
        # BetaD.S[a] * dens.hat[t] +
        # BetaV.S[a] * veg.hat[t] +
        # BetaDV.S[a] * (dens.hat[t] * veg.hat[t]) +
        # BetaVR.S[a] * (veg.hat[t] / dens.hat[t]) +
        Gamma.S[t, a]
    }
  }
  
  # # missing environment
  # for(t in 1:(nYear-1)){
  #   dens.hat[t] ~ dnorm(dens[t], sd = densE[t])
  #   veg.hat[t]  ~ dnorm(veg[t], sd = vegE[t])
  # }
  # 
  # for(m in 1:nNoVeg){
  #   veg[m] ~ dnorm(0, sd = 2)
  # }
  # 
  # for(m in 1:nNoDens){
  #   dens[m] ~ dnorm(0, sd = 2)
  # }
  
  # observation function
  for(t in 1:nYear){
    Epsilon.O[t] ~ dnorm(0, sd = Sigma.O)
    logit(O[t]) <- logit(Mu.O) + Epsilon.O[t]
  }
  
  #### Priors ####
  # for fixed effects
  for(a in 1:nAgeC){
    BetaA.S[a] ~ dunif(-5, 5)
    # BetaD.S[a] ~ dunif(-5, 5)
    # BetaV.S[a] ~ dunif(-5, 5)
    # BetaDV.S[a] ~ dunif(-5, 5)
  }
  
  # for random effects
  # variance-covariance matrix
  for(i in 1:nAgeC){
    zero[i] <- 0
    Xi.S[i] ~ dunif(0, 2)
  }
  
  for(t in 1:(nYear-1)){
    Epsilon.S[t, 1:nAgeC] ~ dmnorm(zero[1:nAgeC], Tau.S[1:nAgeC, 1:nAgeC])
    for(i in 1:nAgeC){
      Gamma.S[t, i] <- Xi.S[i] * Epsilon.S[t, i]
    }
  }
  
  # precision matrix
  Tau.S[1:nAgeC, 1:nAgeC] ~ dwish(W[1:nAgeC, 1:nAgeC], DF)
  Sigma.S[1:nAgeC, 1:nAgeC] <- inverse(Tau.S[1:nAgeC, 1:nAgeC])
  
  # observation
  Mu.O ~ dunif(0.01, 0.99) # or dunif(0, 1)
  Sigma.O ~ dunif(0.01, 10) # or dunif(0, 10)
  
})


## Assemble --------------------------------------------------------------------

nchains   <- 4
seedMod   <- 1:nchains
seedInits <- 1

Tau.S = diag(myConst$nAgeC) + rnorm(myConst$nAgeC^2, 0, 0.1)
Tau.S = inverse((Tau.S + t(Tau.S))/2)

set.seed(seedInits)
myInits <- list()
for(c in 1:nchains){
myInits[[c]] <- list(
  dens = ifelse(is.na(myData$dens), rnorm(myConst$nYear-1, 0, .1), myData$dens),
  veg = ifelse(is.na(myData$veg), rnorm(myConst$nYear-1, 0, .1), myData$veg),
  dens.hat = ifelse(is.na(myData$dens), rnorm(myConst$nYear-1, 0, .1), myData$dens),
  veg.hat = ifelse(is.na(myData$veg), rnorm(myConst$nYear-1, 0, .1), myData$veg),
  
  BetaA.S = rnorm(myConst$nAgeC, 0, 1),
  
  O = runif(myConst$nYear, 0.1, 0.9),
  
  Mu.O = runif(1, 0.1, 0.9),
  Epsilon.O = rnorm(myConst$nYear, 0, 0.2),
  Sigma.O = runif(1, 0.01, 2),
  
  Xi.S = rnorm(myConst$nAgeC, 1, 0.1),
  Epsilon.S = matrix(rnorm((myConst$nYear-1) * myConst$nAgeC, 0, 0.1),
                     nrow = myConst$nYear-1,
                     ncol = myConst$nAgeC),
  
  Tau.S = Tau.S
)}

# select parameters to monitors
params <- c('S', 'BetaA.S', # 'BetaD.S', 'BetaV.S',
            'Mu.O', 'Epsilon.O', 'Sigma.O',
            'Gamma.S', 'Xi.S', 'Sigma.S')

# select MCMC settings
if(testRun){
  nthin   <- 1
  nburnin <- 0
  niter   <- 10
}else{
  nthin   <- 4
  nburnin <- 20000
  niter   <- nburnin + 1000*nthin
}


## Run model -------------------------------------------------------------------

if(parallelRun){
  # function to run one chain inside cluster
  runChain <- function(chainID, code, data, const, inits, params,
                       niter, nburnin, nthin, seed){
    
    library(nimble)
    set.seed(seed)
    inits <- myInits[[chainID]]
    
    model <- nimbleModel(code = myCode, data = myData, constants = myConst, inits = inits)
    cModel <- compileNimble(model)
    conf <- configureMCMC(model, monitors = params)
    mcmc <- buildMCMC(conf)
    cMCMC <- compileNimble(mcmc, project = model)
    
    samples <- runMCMC(cMCMC,
                       thin = nthin,
                       nburnin = nburnin,
                       niter = niter,
                       setSeed = seed,
                       samplesAsCodaMCMC = TRUE)
    return(samples)
  }
  
  # create a cluster & export everything needed to each worker
  cl <- makeCluster(nchains)
  clusterExport(cl, varlist = c("myCode", "myData", "myConst", "myInits", 
                                "params", "nthin", "nburnin", "niter",
                                "seedMod", "runChain"))
}

if(parallelRun){
  # run chains in parallel
  start <- Sys.time()
  samples <- parLapply(cl, 1:nchains, function(i){
    runChain(i,
             code = myCode,
             data = myData,
             const = myConst,
             inits = myInits,
             params = params,
             nthin = nthin,
             nburnin = nburnin, 
             niter = niter, 
             seed = seedMod[i])})
  dur <- Sys.time() - start; dur
  stopCluster(cl)
  beep(2)
}else{
  # run chains sequentially
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
}

# combine & save
out.mcmc <- mcmc.list(samples)
# saveRDS(out.mcmc, 'results/CJS_Age20.rds', compress = 'xz')


## Checks ----------------------------------------------------------------------

library(coda)
library(MCMCvis)
library(corrplot)
library(ggplot2)
library(scales)
library(patchwork)

# # load results
# out.mcmc <- readRDS('results/CJS_Age20.rds')
# summary(out.mcmc) # cannot handle NAs

# summaries
MCMCsummary(out.mcmc, params = c('S', 'BetaA.S'), n.eff = TRUE, round = 2)
# MCMCsummary(out.mcmc, params = c('BetaD.S', 'BetaV.S'), n.eff = TRUE, round = 2)
MCMCsummary(out.mcmc, params = c('Mu.O', 'Epsilon.O', 'Sigma.O'), n.eff = TRUE, round = 2)
MCMCsummary(out.mcmc, params = c('Sigma.S'), n.eff = TRUE, round = 2)

# chainplots
MCMCtrace(out.mcmc, params = c('S', 'BetaA.S'), pdf = FALSE)
# MCMCtrace(out.mcmc, params = c('BetaD.S', 'BetaV.S'), pdf = FALSE)
MCMCtrace(out.mcmc, params = c('Mu.O', 'Epsilon.O', 'Sigma.O'), pdf = FALSE)
MCMCtrace(out.mcmc, params = c('Sigma.S'), pdf = FALSE)


## Plots -----------------------------------------------------------------------

out.dat <- out.mcmc %>% map(as.data.frame) %>% bind_rows()

# check for correlations among fixed effects
par(mfrow = c(1,1))
corrplot(cor(out.mcmc[, grepl('Beta', colnames(out.mcmc[[1]]))] %>%
               map(as.data.frame) %>% bind_rows(), use = 'p'))

# check random effects among demographic rates
# check variance-correlation matrix, with Sigma.S on diagonal
varCorrMatrix <- array(NA, dim = c(myConst$nAgeC, myConst$nAgeC, nrow(out.dat)))

for(i in 1:myConst$nAgeC){
  varCorrMatrix[i,i,] <- out.dat[, paste0('Xi.S[', i,']')]*
    sqrt(out.dat[, paste0('Sigma.S[', i,', ', i,']')])
}

for(j in 1:(myConst$nAgeC-1)){
  for(i in (j+1):myConst$nAgeC){
    varCorrMatrix[j, i, ] <- (out.dat[, paste0('Sigma.S[', i, ', ', j, ']')])/
      sqrt(out.dat[, paste0('Sigma.S[', j, ', ', j, ']')]*
             out.dat[, paste0('Sigma.S[', i, ', ', i, ']')])
  }
}

round(apply(varCorrMatrix, 1:2, mean, na.rm = T), 2)
round(apply(varCorrMatrix, 1:2, quantile, prob = 0.025, na.rm = T), 2)
round(apply(varCorrMatrix, 1:2, quantile, prob = 0.975, na.rm = T), 2)

# calculate survival probabilities
df <- expand.grid(age = 1:20, year = 1:16)
S.pred <- matrix(NA, nrow = nrow(df), ncol = nrow(out.dat))

for(i in 1:nrow(df)){
  S.pred[i, ] <- out.dat[, paste0('BetaA.S[', myData$ageC[df$age[i]], ']')] +
    out.dat[, paste0('Gamma.S[', df$year[i], ', ', myData$ageC[df$age[i]], ']')]
  df$ageC[i] <- myData$ageC[df$age[i]]
}

df$S = inv.logit(apply(S.pred, 1, mean))
df$sLCI = inv.logit(apply(S.pred, 1, quantile, 0.025))
df$sUCI = inv.logit(apply(S.pred, 1, quantile, 0.975))

# plot main results
df %>%
  mutate(ageC = as.factor(ageC)) %>% 
  ggplot(aes(x = year, y = S)) +
  geom_ribbon(aes(ymin = sLCI, ymax = sUCI, fill = ageC), alpha = 0.1) +
  geom_line(aes(colour = ageC), linewidth = 1, show.legend = F) +
  labs(x = "Year", y = "Survival", fill = "Age class") +
  # scale_x_continuous(breaks = pretty_breaks()) +
  scale_y_continuous(breaks = pretty_breaks()) +
  theme_bw()

# ggsave("figures/IPM_CJS.jpeg", scale = 1, width = 18.0, height = 9.0, units = c("cm"), dpi = 600)

young <- df %>%
  mutate(age = as.factor(age-1)) %>% 
  filter(age %in% c(0, 2, 4, 6, 8)) %>%
  ggplot(aes(x = year, y = S)) +
  geom_ribbon(aes(ymin = sLCI, ymax = sUCI, fill = age), alpha = 0.1, show.legend = F) +
  geom_line(aes(colour = age), linewidth = 1) +
  labs(x = "Year", y = "Survival", fill = "Age") +
  # scale_x_continuous(breaks = pretty_breaks()) +
  scale_y_continuous(breaks = pretty_breaks()) +
  theme_bw(); young

old <- df %>%
  mutate(age = as.factor(age-1)) %>% 
  filter(age %in% c(10, 12, 14, 16, 18)) %>%
  ggplot(aes(x = year, y = S)) +
  geom_ribbon(aes(ymin = sLCI, ymax = sUCI, fill = age), alpha = 0.1, show.legend = F) +
  geom_line(aes(colour = age), linewidth = 1) +
  labs(x = "Year", y = "Survival", fill = "Age") +
  # scale_x_continuous(breaks = pretty_breaks()) +
  scale_y_continuous(breaks = pretty_breaks()) +
  theme_bw(); old

young / old

# ggsave("figures/CJS_AgeYr.jpeg", scale = 1, width = 18.0, height = 18.0, units = c("cm"), dpi = 600)

df %>% 
  group_by(age) %>% 
  mutate(age = age-1,
         S = mean(S),
         sUCI = mean(sUCI),
         sLCI = mean(sLCI)) %>% 
  ggplot(aes(x = age, y = S)) +
  geom_ribbon(aes(ymin = sLCI, ymax = sUCI), alpha = 0.1) +
  geom_line(linewidth = 1, show.legend = F) +
  scale_y_continuous(breaks = pretty_breaks()) +
  labs(x = "Age", y = "Survival") +
  theme_bw()

# ggsave("figures/CJS_Age20.jpeg", scale = 1, width = 18.0, height = 9.0, units = c("cm"), dpi = 600)

