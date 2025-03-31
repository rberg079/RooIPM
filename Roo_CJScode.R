# 31 March 2025
# Run survival model

## Set up ----------------------------------------------------------------------

# load packages
library(tidyverse)
library(lubridate)
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

## Model -----------------------------------------------------------------------

myCode = nimbleCode({
  ##### 1. Survival ####
  # Survival function
  for (i in 1:nind){                               
    for (t in first[i]:(last[i]-1)){
      logit(s[i,t]) <- B.age[ageC[age[i,t]]] +
        # veg.hat[t]*B.veg[ageC[age[i,t]]] +
        # dens.hat[t]*B.dens[ageC[age[i,t]]] +
        # (dens.hat[t]*veg.hat[t])*B.densVeg[ageC[age[i,t]]] +
        # (veg.hat[t]/dens.hat[t])*B.vegRoo[ageC[age[i,t]]] +
        gamma[t,ageC[age[i,t]]]
    } #t
  } #i
  
  # Function to estimate missing ages
  for (i in 1:nNoAge){
    ageM[i] ~ T(dnegbin(0.25,1.6),3,20)
    age[noAge[i],first[noAge[i]]] <- round(ageM[i])
    for (t in (first[noAge[i]]+1):ntimes){
      age[noAge[i],t] <- age[noAge[i],t-1]+1
    } #t
  } #i
  
  # Priors for missing values
  # for(t in 1:ntimes){
  #   veg.hat[t] ~ dnorm(veg[t], sd = sd.veg[t])
  #   dens.hat[t] ~ dnorm(dens[t], sd = sd.dens[t])
  # }
  # 
  # for(mt in 1:nNoVeg){
  #   veg[mt] <- 0
  # }
  
  # for(mt in 1:nNoDens){
  #   dens[mt] <- 0
  # }
  
  # Priors for fixed effects
  for(a in 1:nAge){
    B.age[a] ~ dlogis(0,1)
    # B.veg[a] ~ dnorm(0,0.001)
    # B.dens[a] ~ dnorm(0,0.001)
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
  noAge = myconst$noAge
  nNoAge = myconst$nNoAge

  first = myconst$first
  last = myconst$last 

  # veg = mydata$veg
  # dens = mydata$dens
  # nNoVeg = myconst$nNoVeg
  
  # assign initial values
  myinits <- function(i){
    l = list(ageM = sample(3:8, size = nNoAge, replace = T),
             
             B.age = rnorm(nAge,0,0.25),
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
  vars = c('year.p', 'mean.p', 'sd.p', 'state', 'ageM',
           # 'veg.hat', 'sd.veg', 'dens.hat', 'sd.dens',
           'B.age', # 'B.veg', 'B.dens', 'B.densVeg',
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
chain_output <- parLapply(cl = this_cluster,
                          X = 1:3,
                          fun = paraNimble,
                          n.burn = 10000,
                          n.tin = 10,
                          myCode = myCode,
                          myconst = myconst,
                          mydata = mydata)


stopCluster(this_cluster)
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
write_rds(fit1, 'out_fit1.rds', compress = 'xz')


## Checks ----------------------------------------------------------------------

summary(codaSamp)

library(MCMCvis)
MCMCsummary(codaSamp, params = c('B.age','B.veg','B.dens','B.densVeg'), n.eff = TRUE, round = 3)
MCMCsummary(codaSamp, params = c('year.p','mean.p','sd.p'), n.eff = TRUE, round = 3)
MCMCsummary(codaSamp, params = c('Sigma.raw'), n.eff = TRUE, round = 3)
MCMCsummary(codaSamp, params = c('ageM'), n.eff = TRUE, round = 3)

# assess MCMC convergence
# ...visually
par(mar=c(1,1,1,1))
plot(codaSamp[,paste0('B.age[',1:nAge,']')])
plot(codaSamp[,paste0('B.veg[',1:nAge,']')])
plot(codaSamp[,paste0('B.dens[',1:nAge,']')])
plot(codaSamp[,paste0('B.densVeg[',1:nAge,']')])

plot(codaSamp[,paste0('year.p[',1:ntimes,']')])
plot(codaSamp[,'mean.p'])
plot(codaSamp[,'sd.p'])

plot(codaSamp[,paste0('Sigma.raw[',1:nAge,', ',1:nAge,']')])
plot(codaSamp[,paste0('ageM[',1:nb.noAge,']')])

# ...formally
gelman.diag(codaSamp[,paste0('B.age[',1:nAge,']')])
gelman.diag(codaSamp[,paste0('B.veg[',1:nAge,']')])
gelman.diag(codaSamp[,paste0('B.dens[',1:nAge,']')])
gelman.diag(codaSamp[,paste0('B.densVeg[',1:nAge,']')])

gelman.diag(codaSamp[,paste0('year.p[',1:ntimes,']')])
gelman.diag(codaSamp[,'mean.p'])
gelman.diag(codaSamp[,'sd.p'])

gelman.diag(codaSamp[,paste0('Sigma.raw[',1:nAge,', ',1:nAge,']')])
gelman.diag(codaSamp[,paste0('ageM[',1:nb.noAge,']')])

# check Neff
effectiveSize(codaSamp[,paste0('B.age[',1:nAge,']')])
effectiveSize(codaSamp[,paste0('B.veg[',1:nAge,']')])
effectiveSize(codaSamp[,paste0('B.dens[',1:nAge,']')])
effectiveSize(codaSamp[,paste0('B.densVeg[',1:nAge,']')])

effectiveSize(codaSamp[,paste0('year.p[',1:ntimes,']')])
effectiveSize(codaSamp[,'mean.p'])
effectiveSize(codaSamp[,'sd.p'])

effectiveSize(codaSamp[,paste0('Sigma.raw[',1:nAge,', ',1:nAge,']')])
effectiveSize(codaSamp[,paste0('ageM[',1:nb.noAge,']')])


## Plots -----------------------------------------------------------------------

betaz = codaSamp %>% map(as.data.frame) %>% bind_rows()

# check for correlations among fixed effects
tmp <- codaSamp[,grepl('B.',colnames(codaSamp[[1]]))] %>% map(as.data.frame) %>% bind_rows()
corrplot::corrplot(cor(tmp,use = 'p'))

# check random effects among demographic rates
# check variance-correlation matrix, with Sigma.raw on diagonal
VarCorrMatrix <- array(NA,dim = c(myconst$nAge,myconst$nAge,nrow(betaz)))

for (i in 1:myconst$nAge){
  VarCorrMatrix[i,i,] <- betaz[,paste0('xi[',i,']')]*
    sqrt(betaz[,paste0('Sigma.raw[',i,', ',i,']')])
}

for (j in 1:(myconst$nAge-1)){
  for (i in (j+1):myconst$nAge){
    VarCorrMatrix[j,i,] <- ( betaz[,paste0('Sigma.raw[',i,', ',j,']')])/
      sqrt(betaz[,paste0('Sigma.raw[',j,', ',j,']')]*
             betaz[,paste0('Sigma.raw[',i,', ',i,']')])  # correlation coefficients
  }
}

round(apply(VarCorrMatrix, 1:2, mean, na.rm = T), 2)
round(apply(VarCorrMatrix, 1:2, quantile, prob = 0.025, na.rm = T), 2)
round(apply(VarCorrMatrix, 1:2, quantile, prob = 0.975, na.rm = T), 2)

# calculate survival probabilities
df = expand.grid(age = 1:22, year = 1:12)
s.pred = matrix(NA, nrow = nrow(df), ncol = nrow(betaz))

for(i in 1:nrow(df)){
  s.pred[i,] <- betaz[,paste0('B.age[',mydata$ageC[df$age[i]],']')] +
    betaz[,paste0('gamma[',df$year[i],', ', mydata$ageC[df$age[i]],']')]
  df$ageC[i] <- mydata$ageC[df$age[i]]
}

df$s = inv.logit(apply(s.pred,1,mean))
df$s.cil = inv.logit(apply(s.pred,1,quantile,0.025))
df$s.cih = inv.logit(apply(s.pred,1,quantile,0.975))

# plot main results
df <- df %>% mutate(ageC = as.factor(ageC))

library(ggplot2)
library(scales)
plot <- ggplot(df, aes(x = year, y = s)) +
  geom_path(aes(colour = ageC)) +
  labs(x = "Year", y = "Survival", colour = "Age class") +
  scale_x_continuous(breaks = pretty_breaks()) +
  scale_y_continuous(breaks = pretty_breaks(), limits = c(0,1))
plot

