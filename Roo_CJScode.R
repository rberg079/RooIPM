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
source("wrangleData_env.R")
source("wrangleData_surv.R")

env <- wrangleData_env(dens.data = "data/WPNP_Methods_Results_January2025.xlsx",
                       veg.data  = "data/biomass data April 2009 - Jan 2025_updated Feb2025.xlsx",
                       wea.data  = "data/Prom_Weather_2008-2023_updated Jan2025 RB.xlsx",
                       wind.data = "data/POWER_Point_Daily_20080101_20241231_10M.csv")

surv <- wrangleData_surv(surv.data = "data/PromSurvivalOct24.xlsx",
                         yafs.data = "data/RSmainRB_Mar25.xlsx")

# create Nimble lists
myData  <- list(obs = surv$obs,
                state = surv$state,
                age = surv$age,
                ageC = surv$ageC,
                dens = env$dens,
                densE = env$densE,
                veg = env$veg,
                vegE = env$vegE)

myConst <- list(N.id = surv$N.id,
                N.year = surv$N.year,
                N.age = surv$N.age,
                noAge = surv$noAge,
                N.noAge = surv$N.noAge,
                N.noVeg = env$N.noVeg,
                first = surv$first,
                last = surv$last,
                W = surv$W,
                DF = surv$DF)

# # checks
# sapply(1:nrow(state), function (i) state[i, first[i]]) %>% table(useNA = 'a') # should all be 1
# sapply(1:nrow(state), function (i) state[i, last[i]]) %>% table(useNA = 'a')  # should be mostly 0s & NAs
# sapply(1:nrow(state), function (i) age[i, first[i]]) %>% table(useNA = 'a')   # should all be >= 1
# table(first >= last, useNA = 'a')                                             # should all be F

## Model -----------------------------------------------------------------------

myCode = nimbleCode({
  ##### 1. Survival ####
  # Survival function
  for (i in 1:N.id){                               
    for (t in first[i]:(last[i]-1)){
      logit(s[i, t]) <- B.age[ageC[age[i, t]]] +
        B.dens[ageC[age[i, t]]] * dens.hat[t] +
        B.veg[ageC[age[i, t]]] * veg.hat[t] +
        # B.densVeg[ageC[age[i, t]]] * (dens.hat[t] * veg.hat[t]) +
        # B.vegRoo[ageC[age[i, t]]] * (veg.hat[t] / dens.hat[t]) +
        gamma[t, ageC[age[i,t]]]
    } # t
  } # i
  
  for (t in 1:N.year){
    dens.hat[t] ~ dnorm(dens[t], sd = densE[t])
    veg.hat[t] ~ dnorm(veg[t], sd = vegE[t])
  }
  
  for (mt in 1:N.noVeg){
    veg[mt] ~ dnorm(0,0.001)
  }
  
  # Estimate missing ages
  for (i in 1:N.noAge){
    ageM[i] ~ T(dnegbin(0.25, 1.6), 3, 20)
    age[noAge[i], first[noAge[i]]] <- round(ageM[i]) + 1
    for (t in (first[noAge[i]]+1):N.year){
      age[noAge[i], t] <- age[noAge[i], t-1] + 1
    } # t
  } # i
  
  # Priors for fixed effects
  for(a in 1:N.age){
    B.age[a] ~ dlogis(0, 1)
    B.dens[a] ~ dnorm(0, 0.001)
    B.veg[a] ~ dnorm(0, 0.001)
    # B.densVeg[a] ~ dnorm(0, 0.001)
  }
  
  # Variance-Covariance matrix
  for (i in 1:N.age){
    zero[i] <- 0
    xi[i] ~ dunif(0, 2)
  } # i
  
  for (t in 1:(N.year-1)){
    eps.raw[t, 1:N.age]  ~ dmnorm(zero[1:N.age], Tau.raw[1:N.age, 1:N.age])
    for (i in 1:N.age){
      gamma[t, i] <- xi[i] * eps.raw[t,i]
    } # i
  } # t
  
  # Priors for precision matrix
  Tau.raw[1:N.age, 1:N.age] ~ dwish(W[1:N.age, 1:N.age], DF)
  Sigma.raw[1:N.age, 1:N.age] <- inverse(Tau.raw[1:N.age, 1:N.age])
  
  # Uniform covariance matrix
  # for (a in 1:N.age){
  #   zero[a] <- 0
  #   sd.yr[a] ~ dunif(0, 5)
  #   cov.yr[a, a] <- sd.yr[a] * sd.yr[a]
  # }
  #
  # for(a in 1:(N.age-1)){
  #   for(a2 in (a+1):N.age){
  #     cor.yr[a, a2] ~ dunif(-1, 1)
  #     cov.yr[a2, a] <- sd.yr[a] * sd.yr[a2] * cor.yr[a, a2]
  #     cov.yr[a, a2] <- cov.yr[a2, a]
  #   }
  # } # i
  #
  # for (t in 1:(N.year-1)) {
  #   gamma[t, 1:N.age]  ~ dmnorm(zero[1:N.age], cov = cov.yr[1:N.age, 1:N.age])
  # } # t
  
  
  ##### 2. Observation ####
  for (t in 1:N.year){
    for(i in 1:N.id){
      logit(p[i, t]) <- mu.p + year.p[t]
    }
    year.p[t] ~ dnorm(0, sd = sd.p)
  }
  
  mu.p <- log(mean.p / (1-mean.p))
  mean.p ~ dunif(0, 1)
  sd.p ~ dunif(0, 10) # try dunif(0.01, 10) if sd.p causes trouble
  
  
  ##### 3. Likelihood ####
  for (i in 1:N.id){
    for (t in (first[i] + 1):last[i]){
      # State process
      state[i, t] ~ dbern(mu1[i, t])
      mu1[i, t] <- s[i, t-1] * state[i, t-1]
      
      # Observation process
      obs[i, t] ~ dbern(mu2[i, t])
      mu2[i, t] <- p[i, t] * state[i, t]
    } #t
  } #i
  
})


## Assemble --------------------------------------------------------------------

# create Nimble function
paraNimble <- function(seed, myCode, myConst, myData,
                       n.burn = 1000, n.tin = 1){
  
  library(nimble)
  library(coda)
  
  N.id = myConst$N.id
  N.year = myConst$N.year
  
  age = myData$age
  N.age = myConst$N.age
  noAge = myConst$noAge
  N.noAge = myConst$N.noAge
  
  first = myConst$first
  last = myConst$last 
  
  veg = myData$veg
  dens = myData$dens
  N.noVeg = myConst$N.noVeg
  
  # assign initial values
  myInits <- function(i){
    l = list(ageM = sample(3:8, size = N.noAge, replace = T),
             
             B.age = rnorm(N.age, 0, 0.25),
             B.veg = rnorm(N.age, 0, 0.5),
             B.dens = rnorm(N.age, 0, 1),
             # B.densVeg = rnorm(N.age, 0, 1),
             
             veg.hat = ifelse(is.na(veg), rnorm(length(veg), 0, .1), veg),
             dens.hat = ifelse(is.na(dens), rnorm(length(dens), 0, .1), dens),
             
             mean.p = runif(1, 0.6, 1),
             year.p = rnorm(N.year, 0, 0.2),
             sd.p = rnorm(1, 0.2, 0.1),
             
             xi = rnorm(N.age, 1, 0.1),
             eps.raw = matrix(rnorm((N.year-1)*N.age, 0, 0.1),
                              ncol = N.age, nrow = (N.year-1))
             
             # cor.yr = diag(N.age) + 0.01,
             # sd.yr = runif(N.age, 0,1)
    )
    Tau.raw = diag(N.age) + rnorm(N.age^2, 0, 0.1)
    l$Tau.raw = inverse((Tau.raw + t(Tau.raw))/2)  # may be a typo?
    return(l)
  }
  
  # assemble model
  myMod <- nimbleModel(code = myCode,
                       data = myData,
                       constants = myConst,
                       inits = myInits())
  
  # select parameters to monitor
  vars = c('year.p', 'mean.p', 'sd.p', 'state', 'ageM',
           'B.age', 'B.dens', 'B.veg', # 'B.densVeg',
           'dens.hat', 'veg.hat', 'densE', 'vegE',
           'gamma', 'xi', 'Sigma.raw'
           # 'gamma', 'sd.yr', 'cor.yr'
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

start.t <- Sys.time()
this_cluster <- makeCluster(3)
chain_output <- parLapply(X = 1:3,
                          cl = this_cluster,
                          fun = paraNimble,
                          myCode = myCode,
                          myconst = myConst,
                          mydata = myData,
                          n.burn = 10000,
                          n.tin = 10)

beep(sound = 2)
stopCluster(this_cluster)
dur = now() - start.t
dur

# reformat output
codaSamp <- chain_output %>% map(~ as.mcmc(.x$samples)) %>% as.mcmc.list()
codaSamp <- codaSamp[, !grepl('state', colnames(codaSamp[[1]]))]

# obtain WAIC value
waic <- map_dbl(chain_output, ~ .x$WAIC$WAIC)  # unserialized
# waic <- chain_output %>% map_dbl(~.x$WAIC)   # serialized
waic

# save output
fit1 <- list(model = myCode, codaSamp = codaSamp, waic = waic, dur = dur)
# write_rds(fit1, 'results/out_fit1.rds', compress = 'xz')


## Checks ----------------------------------------------------------------------

summary(codaSamp) # cannot handle any NAs

library(MCMCvis)
library(corrplot)
MCMCsummary(codaSamp, params = c('B.age', 'B.dens', 'B.veg'), n.eff = TRUE, round = 3) # 'B.densVeg'
MCMCsummary(codaSamp, params = c('year.p','mean.p','sd.p'), n.eff = TRUE, round = 3)
MCMCsummary(codaSamp, params = c('Sigma.raw'), n.eff = TRUE, round = 3)
MCMCsummary(codaSamp, params = c('ageM'), n.eff = TRUE, round = 3)

# assess MCMC convergence
# ...visually
par(mar = c(1,1,1,1))
plot(codaSamp[, paste0('B.age[',1:N.age,']')])
plot(codaSamp[, paste0('B.veg[',1:N.age,']')])
plot(codaSamp[, paste0('B.dens[',1:N.age,']')])
plot(codaSamp[, paste0('B.densVeg[',1:N.age,']')])

plot(codaSamp[, paste0('year.p[',1:N.year,']')])
plot(codaSamp[, 'mean.p'])
plot(codaSamp[, 'sd.p'])

plot(codaSamp[, paste0('Sigma.raw[',1:N.age,', ',1:N.age,']')])
plot(codaSamp[, paste0('ageM[',1:nb.noAge,']')])

# # ...formally
# gelman.diag(codaSamp[, paste0('B.age[',1:N.age,']')])
# gelman.diag(codaSamp[, paste0('B.veg[',1:N.age,']')])
# gelman.diag(codaSamp[, paste0('B.dens[',1:N.age,']')])
# gelman.diag(codaSamp[, paste0('B.densVeg[',1:N.age,']')])
# 
# gelman.diag(codaSamp[, paste0('year.p[',1:N.year,']')])
# gelman.diag(codaSamp[, 'mean.p'])
# gelman.diag(codaSamp[, 'sd.p'])
# 
# gelman.diag(codaSamp[, paste0('Sigma.raw[',1:N.age,', ',1:N.age,']')])
# gelman.diag(codaSamp[, paste0('ageM[',1:nb.noAge,']')])
# 
# # check Neff
# effectiveSize(codaSamp[, paste0('B.age[',1:N.age,']')])
# effectiveSize(codaSamp[, paste0('B.veg[',1:N.age,']')])
# effectiveSize(codaSamp[, paste0('B.dens[',1:N.age,']')])
# effectiveSize(codaSamp[, paste0('B.densVeg[',1:N.age,']')])
# 
# effectiveSize(codaSamp[, paste0('year.p[',1:N.year,']')])
# effectiveSize(codaSamp[, 'mean.p'])
# effectiveSize(codaSamp[, 'sd.p'])
# 
# effectiveSize(codaSamp[, paste0('Sigma.raw[',1:N.age,', ',1:N.age,']')])
# effectiveSize(codaSamp[, paste0('ageM[',1:nb.noAge,']')])


## Plots -----------------------------------------------------------------------

betaz = codaSamp %>% map(as.data.frame) %>% bind_rows()

# check for correlations among fixed effects
par(mfrow = c(1,1))
corrplot(cor(codaSamp[, grepl('B.', colnames(codaSamp[[1]]))] %>%
               map(as.data.frame) %>%
               bind_rows()
             , use = 'p'))

# check random effects among demographic rates
# check variance-correlation matrix, with Sigma.raw on diagonal
varCorrMatrix <- array(NA, dim = c(myconst$N.age, myconst$N.age, nrow(betaz)))

for (i in 1:myconst$N.age){
  varCorrMatrix[i,i,] <- betaz[, paste0('xi[', i,']')]*
    sqrt(betaz[, paste0('Sigma.raw[', i,', ', i,']')])
}

for (j in 1:(myconst$N.age-1)){
  for (i in (j+1):myconst$N.age){
    varCorrMatrix[j,i,] <- (betaz[, paste0('Sigma.raw[',i,', ',j,']')])/
      sqrt(betaz[, paste0('Sigma.raw[', j,', ', j,']')]*
             betaz[, paste0('Sigma.raw[', i,', ', i,']')])
  }
}

round(apply(varCorrMatrix, 1:2, mean, na.rm = T), 2)
round(apply(varCorrMatrix, 1:2, quantile, prob = 0.025, na.rm = T), 2)
round(apply(varCorrMatrix, 1:2, quantile, prob = 0.975, na.rm = T), 2)

# calculate survival probabilities
df = expand.grid(age = 1:22, year = 1:16)
s.pred = matrix(NA, nrow = nrow(df), ncol = nrow(betaz))

for(i in 1:nrow(df)){
  s.pred[i,] <- betaz[, paste0('B.age[', mydata$ageC[df$age[i]], ']')] +
    betaz[,paste0('gamma[', df$year[i], ', ', mydata$ageC[df$age[i]], ']')]
  df$ageC[i] <- mydata$ageC[df$age[i]]
}

df$s = inv.logit(apply(s.pred, 1, mean))
df$s.cil = inv.logit(apply(s.pred, 1, quantile, 0.025))
df$s.cih = inv.logit(apply(s.pred, 1, quantile, 0.975))

# plot main results
df <- df %>% mutate(ageC = as.factor(ageC))

library(ggplot2)
library(scales)
plot <- ggplot(df, aes(x = year, y = s)) +
  geom_path(aes(colour = ageC)) +
  labs(x = "Year", y = "Survival", colour = "Age class") +
  # scale_x_continuous(breaks = pretty_breaks()) +
  scale_y_continuous(breaks = pretty_breaks(),
                     limits = c(0,1)) +
  theme_bw(); plot

