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
enData <- wrangleData_env(dens.data = "data/WPNP_Methods_Results_January2025.xlsx",
                          veg.data  = "data/biomass data April 2009 - Jan 2025_updated Feb2025.xlsx",
                          wea.data  = "data/Prom_Weather_2008-2023_updated Jan2025 RB.xlsx",
                          wind.data = "data/POWER_Point_Daily_20080101_20241231_10M.csv")

source("wrangleData_sv.R")
svData <- wrangleData_surv(surv.data = "data/PromSurvivalOct24.xlsx",
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

myConst <- list(nID.sv = svData$nID.sv,
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

## Model -----------------------------------------------------------------------

myCode = nimbleCode({
  
  #### Likelihood ####
  for (i in 1:nID.sv){
    for (t in (first[i] + 1):last[i]){
      # state process
      state[i, t] ~ dbern(MuS[i, t])
      MuS[i, t] <- sv[ageC[age[i, t-1]], t-1] * state[i, t-1]
      # observation process
      obs[i, t] ~ dbern(MuO[i, t])
      MuO[i, t] <- ob[i, t] * state[i, t]
    }
  }
  
  #### Constraints ####
  # survival function
  for (i in 1:nID.sv){                               
    for (t in first[i]:(last[i]-1)){
      logit(sv[i, t]) <- BetaA.sv[ageC[age[i, t]]] +
        BetaD.sv[ageC[age[i, t]]] * dens.hat[t] +
        BetaV.sv[ageC[age[i, t]]] * veg.hat[t] +
        # BetaDV.sv[ageC[age[i, t]]] * (dens.hat[t] * veg.hat[t]) +
        # BetaVR.sv[ageC[age[i, t]]] * (veg.hat[t] / dens.hat[t]) +
        Gamma.sv[t, ageC[age[i,t]]]
    }
  }
  
  for (t in 1:nYear){
    dens.hat[t] ~ dnorm(dens[t], sd = densE[t])
    veg.hat[t] ~ dnorm(veg[t], sd = vegE[t])
  }
  
  for (m in 1:nNoVeg){
    veg[m] <- 0
    # veg[m] ~ dnorm(0, 0.001)
  }
  
  # estimate missing ages
  for (i in 1:nNoAge){
    ageM[i] ~ T(dnegbin(0.25, 1.6), 3, 20)
    age[noAge[i], first[noAge[i]]] <- round(ageM[i]) + 1
    for (t in (first[noAge[i]]+1):nYear){
      age[noAge[i], t] <- age[noAge[i], t-1] + 1
    }
  }
  
  # observation function
  for (t in 1:nYear){
    for(i in 1:nID.sv){
      logit(ob[i, t]) <- logit(MuO) + Epsilon.ob[t]
    }
    Epsilon.ob[t] ~ dnorm(0, sd = Sigma.ob)
  }
  
  #### Priors ####
  # for fixed effects
  for(a in 1:nAgeC){
    BetaA.sv[a] ~ dlogis(0, 1) # TODO: think about this?
    BetaD.sv[a] ~ dnorm(0, 0.001)
    BetaV.sv[a] ~ dnorm(0, 0.001)
    # BetaDV.sv[a] ~ dnorm(0, 0.001)
  }
  
  # for random effects
  # variance-covariance matrix
  for (i in 1:nAgeC){
    zero[i] <- 0
    Xi.sv[i] ~ dunif(0, 2)
  }
  
  for (t in 1:(nYear-1)){
    Epsilon.sv[t, 1:nAgeC] ~ dmnorm(zero[1:nAgeC], Tau.sv[1:nAgeC, 1:nAgeC])
    for (i in 1:nAgeC){
      Gamma.sv[t, i] <- Xi.sv[i] * Epsilon.sv[t,i]
    }
  }
  
  # precision matrix
  Tau.sv[1:nAgeC, 1:nAgeC] ~ dwish(W[1:nAgeC, 1:nAgeC], DF)
  Sigma.sv[1:nAgeC, 1:nAgeC] <- inverse(Tau.sv[1:nAgeC, 1:nAgeC])
  
  # observation
  MuO ~ dunif(0.2, 1)
  Sigma.ob ~ dunif(0.01, 10) # or dunif(0, 10)
  
})


## Assemble --------------------------------------------------------------------

# create Nimble function
paraNimble <- function(seed, myCode, myConst, myData,
                       n.burn = 1000, n.thin = 1){
  
  library(nimble)
  library(coda)
  
  nID.sv = myConst$nID.sv
  nYear = myConst$nYear
  
  age = myData$age
  nAgeC = myConst$nAgeC
  noAge = myConst$noAge
  nNoAge = myConst$nNoAge
  
  first = myConst$first
  last = myConst$last 
  
  veg = myData$veg
  dens = myData$dens
  nNoVeg = myConst$nNoVeg
  
  # assign initial values
  myInits <- function(i){
    l = list(ageM = sample(3:8, size = nNoAge, replace = T),
             
             BetaA.sv = rnorm(nAgeC, 0, 1),
             BetaV.sv = rnorm(nAgeC, 0, 1),
             BetaD.sv = rnorm(nAgeC, 0, 1),
             # BetaDV.sv = rnorm(nAgeC, 0, 1),
             
             veg.hat = ifelse(is.na(veg), rnorm(length(veg), 0, .1), veg),
             dens.hat = ifelse(is.na(dens), rnorm(length(dens), 0, .1), dens),
             
             MuO = runif(1, 0.2, 0.8),
             Epsilon.ob = rnorm(nYear, 0, 0.2),
             Sigma.ob = runif(1, 0.01, 2), # or rnorm(1, 0.2, 0.1)
             
             Xi.sv = rnorm(nAgeC, 1, 0.1),
             Epsilon.sv = matrix(rnorm((nYear-1)*nAgeC, 0, 0.1),
                                 ncol = nAgeC, nrow = (nYear-1))
    )
    Tau.sv = diag(nAgeC) + rnorm(nAgeC^2, 0, 0.1)
    l$Tau.sv = inverse((Tau.sv + t(Tau.sv))/2) # may be a typo?
    return(l)
  }
  
  # assemble model
  myMod <- nimbleModel(code = myCode,
                       data = myData,
                       constants = myConst,
                       inits = myInits())
  
  # select parameters to monitor
  vars = c('Epsilon.ob', 'Sigma.ob', 'state', 'ageM',
           'BetaA.sv', 'BetaD.sv', 'BetaV.sv', # 'BetaDV.sv',
           'dens.hat', 'veg.hat', 'densE', 'vegE',
           'Gamma.sv', 'Xi.sv', 'Sigma.sv'
  )
  
  # select MCMC settings
  cModel <- compileNimble(myMod)
  mymcmc <- buildMCMC(cModel, monitors = vars, enableWAIC = T)
  CmyMCMC <- compileNimble(mymcmc, project = myMod)
  samples <- runMCMC(CmyMCMC,
                     samplesAsCodaMCMC = T,
                     niter = n.burn + 1000*n.thin,
                     nburnin = n.burn,
                     thin = n.thin,
                     summary = T,
                     WAIC = T)
  
  return(samples)
}


## Run model -------------------------------------------------------------------

start.t <- Sys.time()
this_cluster <- makeCluster(3)
samples <- parLapply(X = 1:3,
                     cl = this_cluster,
                     fun = paraNimble,
                     myCode = myCode,
                     myConst = myConst,
                     myData = myData,
                     n.burn = 10000,
                     n.thin = 10)

beep(sound = 2)
stopCluster(this_cluster)
dur = now() - start.t
dur

# reformat output
out.mcmc <- samples %>% map(~ as.mcmc(.x$samples)) %>% as.mcmc.list()
out.mcmc <- out.mcmc[, !grepl('state', colnames(out.mcmc[[1]]))]

# obtain WAIC value
waic <- samples %>% map_dbl(~.x$WAIC)  # serialized
waic

# save output
fit <- list(model = myCode, out.mcmc = out.mcmc, waic = waic, dur = dur)
# write_rds(fit, 'results/out_fit.rds', compress = 'xz')


## Checks ----------------------------------------------------------------------

summary(out.mcmc) # cannot handle any NAs

library(MCMCvis)
library(corrplot)
MCMCsummary(out.mcmc, params = c('BetaA.sv', 'BetaD.sv', 'BetaV.sv'), n.eff = TRUE, round = 3)
MCMCsummary(out.mcmc, params = c('Epsilon.ob','Sigma.ob'), n.eff = TRUE, round = 3)
MCMCsummary(out.mcmc, params = c('Sigma.sv'), n.eff = TRUE, round = 3)
MCMCsummary(out.mcmc, params = c('ageM'), n.eff = TRUE, round = 3)

# assess MCMC convergence
# ...visually
par(mar = c(1,1,1,1))
plot(out.mcmc[, paste0('BetaA.sv[',1:nAgeC,']')])
plot(out.mcmc[, paste0('BetaV.sv[',1:nAgeC,']')])
plot(out.mcmc[, paste0('BetaD.sv[',1:nAgeC,']')])
plot(out.mcmc[, paste0('BetaDV.sv[',1:nAgeC,']')])

plot(out.mcmc[, paste0('Epsilon.ob[',1:nYear,']')])
plot(out.mcmc[, 'Sigma.ob'])

plot(out.mcmc[, paste0('Sigma.sv[',1:nAgeC,', ',1:nAgeC,']')])
plot(out.mcmc[, paste0('ageM[',1:nNoAge,']')])

# # ...formally
# gelman.diag(out.mcmc[, paste0('BetaA.sv[',1:nAgeC,']')])
# gelman.diag(out.mcmc[, paste0('BetaV.sv[',1:nAgeC,']')])
# gelman.diag(out.mcmc[, paste0('BetaD.sv[',1:nAgeC,']')])
# gelman.diag(out.mcmc[, paste0('BetaDV.sv[',1:nAgeC,']')])
# 
# gelman.diag(out.mcmc[, paste0('Epsilon.ob[',1:nYear,']')])
# gelman.diag(out.mcmc[, 'Sigma.ob'])
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
# effectiveSize(out.mcmc[, paste0('Epsilon.ob[',1:nYear,']')])
# effectiveSize(out.mcmc[, 'Sigma.ob'])
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

for (i in 1:myConst$nAgeC){
  varCorrMatrix[i,i,] <- betaz[, paste0('Xi.sv[', i,']')]*
    sqrt(betaz[, paste0('Sigma.sv[', i,', ', i,']')])
}

for (j in 1:(myConst$nAgeC-1)){
  for (i in (j+1):myConst$nAgeC){
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
  geom_path(aes(colour = ageC)) +
  labs(x = "Year", y = "Survival", colour = "Age class") +
  # scale_x_continuous(breaks = pretty_breaks()) +
  scale_y_continuous(breaks = pretty_breaks(),
                     limits = c(0,1)) +
  theme_bw(); plot

