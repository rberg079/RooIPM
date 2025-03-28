
# Supporting information for "Post-weaning survival in kangaroos is high and constant until senescence: implications for population dynamics"
# Rachel Bergeron, Gabriel Pigeon, David M. Forsyth, Wendy J. King, and Marco Festa-Bianchet 2022
# Ecology

# Script to assemble Nimble data and run state-space formulation of Cormack-Jolly-Seber model
# to estimate apparent survival while accounting for imperfect detection

########################################

## Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## Use packages
library(tidyverse)
library(lubridate)
library(boot)
library(coda)
library(foreach)
library(doParallel)
library(parallel)
library(nimble)
registerDoParallel(3)

## Import data
# For female models...
obs <- read_csv('Data/obsF.csv') %>% as.matrix()
state <- read_csv('Data/stateF.csv') %>% as.matrix()
age <- read_csv('Data/ageF.csv') %>% as.matrix()
id <- read_csv('Data/idF.csv')

# ... or for male models
# obs <- read_csv('Data/obsM.csv') %>% as.matrix()
# state <- read_csv('Data/stateM.csv') %>% as.matrix()
# age <- read_csv('Data/ageM.csv') %>% as.matrix()
# id <- read_csv('Data/idM.csv')

# with current...
env <- read_csv('Data/env.csv')

# ... or time-lagged environmental variables
# env <- read_csv('Data/envLag.csv')

# Create age classes
SurvAgeClasses = c(0,1,1,rep(2,4), rep(3,3), rep(4,60))

## Collate data for Nimble
# Individual data
nind = nrow(state)
ntimes = ncol(state)

nAge = max(SurvAgeClasses, na.rm = T)+1
noAge = which(is.na(age[,ncol(age)]))    # Identify individuals of unknown age
nb.noAge = length(noAge)                 # Number of individuals of unknwon age

first = as.numeric(id$first)
last = as.numeric(id$last)

SageC = SurvAgeClasses+1    # +1 so age and age classes start at 1 rather than 0
age = as.matrix(age)+1

# Environmental data

# Because Veg and Dens are highly correlated, when both included
# in the same model, the common variation was subtracted from Dens using
# sequential regression (Dormann et al. 2013) and resDens was used instead of Dens

veg = round(as.numeric(scale(env$Veg)),3)         # Scaled
# dens = round(as.numeric(scale(env$Dens)),3)     # Scaled
dens = round(as.numeric(env$resDens),3)

sd.veg = round(as.numeric(env$VegSD/sd(env$Veg, na.rm = T)),3)
# sd.dens = round(as.numeric(env$DensSD/sd(env$Dens, na.rm = T)),3)
sd.dens = round(as.numeric(env$resDensSD/sd(env$resDens, na.rm = T)),3)

sd.veg = ifelse(is.na(sd.veg), 2, sd.veg)         # Set high SD for years of
sd.dens = ifelse(is.na(sd.dens), 2, sd.dens)      # missing Veg or Dens data

nmissingVeg = sum(is.na(veg))                     # Years of missing Veg data
nmissingDens = sum(is.na(dens))                   # Years of missing Dens data

## Assemble Nimble lists
mydata <- list(age = age, SageC = SageC, obs = obs, state = state,
               veg = veg, sd.veg = sd.veg, dens = dens, sd.dens = sd.dens)

myconst <- list(nind = nind, ntimes = ntimes,
                nAge = nAge, noAge = noAge, nb.noAge = nb.noAge,
                first = first, last = last, W = diag(nAge), DF = nAge+1,
                nmissingVeg = nmissingVeg, nmissingDens = nmissingDens)

## Write Nimble code
myCode = nimbleCode({
  # --------------------------------------
  # 1.Survival
  # --------------------------------------
  # Estimate age at first capture of individuals of unknown age
  for (i in 1:nb.noAge){
    ageM[i] ~ T(dnegbin(0.25,1.6),3,20)                  # Priors drawn from a negative binomial distribution
    age[noAge[i],first[noAge[i]]] <- round(ageM[i])      # approximating the average age distribution of the
    for (t in (first[noAge[i]]+1):ntimes){               # population, truncated because individuals caught
      age[noAge[i],t] <- age[noAge[i],t-1]+1             # very young would be of known age
    } #t
  } #i
  
  # Survival function
  for (i in 1:nind){                               
    for (t in first[i]:(last[i]-1)){
      logit(s[i,t]) <- B.age[SageC[age[i,t]]] +          # Mixed generalized regression with logit link with no intercept
        veg.hat[t]*B.veg[SageC[age[i,t]]] +              # Full model with age-specific effects of Veg, Dens, and Veg*Dens
        dens.hat[t]*B.dens[SageC[age[i,t]]] +
        (dens.hat[t]*veg.hat[t])*B.densVeg[SageC[age[i,t]]] +
        # (veg.hat[t]/dens.hat[t])*B.vegRoo[SageC[age[i,t]]] +
        gamma[t,SageC[age[i,t]]]
    } #t
  } #i
  
  for(t in 1:ntimes){
    veg.hat[t] ~ dnorm(veg[t], sd = sd.veg[t])           # True Veg and Dens were drawn from normal prior
    dens.hat[t] ~ dnorm(dens[t], sd = sd.dens[t])        # distributions centered around their estimates
  }
  
  # Set prior value of 0 (average) for years of missing Veg and Dens data
  for(mt in 1:nmissingVeg){
    veg[mt] <- 0
  }
  
  for(mt in 1:nmissingDens){
    dens[mt] <- 0
  }
  
  # Variance-Covariance matrix
  # Priors for random effects of demographic rates drawn from
  # the scaled inverse Wishart distribution (Gelman and Hill 2007, p 376)
  for (i in 1:nAge){
    zero[i] <- 0
    xi[i] ~ dunif(0,2) # Scaling
  } #i
  for (t in 1:(ntimes-1)){
    eps.raw[t,1:nAge]  ~ dmnorm(zero[1:nAge], Tau.raw[1:nAge,1:nAge])
    for (i in 1:nAge){
      gamma[t,i] <- xi[i] * eps.raw[t,i]
    } #i
  } #t

  # Priors for precision matrix
  Tau.raw[1:nAge, 1:nAge] ~ dwish(W[1:nAge,1:nAge], DF)           # W = diag(nAge) provided in myconst
  Sigma.raw[1:nAge, 1:nAge] <- inverse(Tau.raw[1:nAge,1:nAge])    # DF = nAge+1 provided in myconst
  
  # Uniform covariance matrix
  # Used in a sensitivity analysis to compare to
  # results obtained using the scaled inverse Wishart distribution
  # for (a in 1:nAge){
  #   zero[a] <- 0
  #   sd.yr[a] ~ dunif(0,5) # Scaling
  #   cov.yr[a,a] <- sd.yr[a]*sd.yr[a]
  # }
  # for(a in 1:(nAge-1)){
  #   for(a2 in (a+1):nAge){
  #     cor.yr[a,a2] ~ dunif(-1,1)
  #     cov.yr[a2,a] <- sd.yr[a] * sd.yr[a2] * cor.yr[a,a2]
  #     cov.yr[a,a2] <- cov.yr[a2,a]
  #   }
  # } #i
  # for (t in 1:(ntimes-1)) {
  #   gamma[t,1:nAge]  ~ dmnorm(zero[1:nAge], cov=cov.yr[1:nAge, 1:nAge])
  # } #t
  
  # Priors for fixed effects
  for(a in 1:nAge){
    B.age[a] ~ dlogis(0,1)
    B.veg[a] ~ dnorm(0,0.001)
    B.dens[a] ~ dnorm(0,0.001)
    B.densVeg[a] ~ dnorm(0,0.001)
  }
  
  # --------------------------------------
  # 2.Observation
  # --------------------------------------
  for (t in 1:ntimes){
    for(i in 1:nind){
      logit(p[i,t]) <- mu.p + year.p[t]    # Probability of observing i at time t depends on a mean
    }                                      # observation probability mu.p and yearly variation year.p
    year.p[t] ~ dnorm(0, sd = sd.p)
  }
  
  mu.p <- log(mean.p / (1-mean.p))
  mean.p ~ dunif(0,1)
  sd.p ~ dunif(0,10)
  
  # --------------------------------------
  # 3.Likelihood
  # --------------------------------------
  for (i in 1:nind){
    for (t in (first[i] + 1):last[i]){
      # State process                      # State taken from a Bernoulli distribution
      state[i,t] ~ dbern(mu1[i,t])         # Probability of being alive at time t (mu1) depends on
      mu1[i,t] <- s[i,t-1] * state[i,t-1]  # state at time t-1 and survival probability (s)
  
      # Observation process                # Observation taken from a Bernoulli distribution
      obs[i,t] ~ dbern(mu2[i,t])           # Probability of being observed at time t (mu2) depends on
      mu2[i,t] <- p[i,t] * state[i,t]      # state at time t and observation probability (p)
    } #t
  } #i
  
})

## Create Nimble function for parallel processing
paraNimble <- function(seed, myCode, myconst, mydata,
                       n.burn = 1000, n.tin = 1){
  
  library(nimble)
  library(coda)
  
  nind = myconst$ntimes
  ntimes = myconst$ntimes
  
  nAge = myconst$nAge
  noAge = myconst$noAge
  nb.noAge = myconst$nb.noAge

  first = myconst$first
  last = myconst$last  
  age = mydata$age

  veg = mydata$veg
  dens = mydata$dens
  nmissingVeg = myconst$nmissingVeg
  nmissingDens = myconst$nmissingDens
  
  # Assign initial values
  # Taken randomly around plausible values
  myinits <- function(i){
    l = list(ageM = sample(3:8, size = nb.noAge, replace = T),    # Prior for missing ages
             
             B.age = rnorm(nAge,0,0.25),                          # Priors for fixed effects
             B.veg = rnorm(nAge,0,0.5),
             B.dens = rnorm(nAge,0,1),
             B.densVeg = rnorm(nAge,0,1),
             
             veg.hat = ifelse(is.na(veg),rnorm(length(veg),0,.1),veg),        # Priors for true values
             dens.hat = ifelse(is.na(dens),rnorm(length(dens),0,.1),dens),    # of Veg and Dens
             
             mean.p = runif(1,0.6,1),                             # Priors for mean observation
             year.p = rnorm(ntimes,0,0.2),                        # probability and year effect
             sd.p = rnorm(1,0.2,0.1),
             
             xi = rnorm(nAge,1,0.1),                              # Priors for elements of the 
             eps.raw = matrix(rnorm((ntimes-1)*nAge,0,0.1),       # variance-covariance matrix
                              ncol = nAge, nrow = (ntimes-1))
             
             # cor.yr = diag(nAge)+0.01,                          # Priors for elements of the
             # sd.yr = runif(nAge,0,1)                            # uniform covariance matrix
    )
    Tau.raw = diag(nAge) + rnorm(nAge^2,0,0.1)
    l$Tau.raw = inverse((Tau.raw + t(Tau.raw))/2)
    return(l)
  }
  
  # Assemble model
  myMod <- nimbleModel(code = myCode,
                       data = mydata,
                       constants = myconst,
                       inits = myinits())
  
  # Select variables to be monitored
  vars = c('year.p', 'mean.p', 'sd.p', 'state', 'ageM',
           'veg.hat', 'sd.veg', 'dens.hat', 'sd.dens',
           'B.age', 'B.veg', 'B.dens', 'B.densVeg',
           'gamma', 'xi', 'Sigma.raw'                      # When using variance-covariance matrix
           # 'gamma', 'sd.yr', 'cor.yr'                    # When using uniform covariance matrix
  )
  
  # Select MCMC settings
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

## Run the model
start.t <- Sys.time()
this_cluster <- makeCluster(3)
chain_output <- parLapply(cl = this_cluster, X = 1:3,      # Run 3 chains
                          fun = paraNimble,
                          n.burn = 10000,         # Increase for serious inference
                          n.tin = 10,
                          myCode = myCode,
                          myconst = myconst,
                          mydata = mydata)


stopCluster(this_cluster)
dur = now() - start.t
dur

## Reformat output
codaSamp <- chain_output %>% map(~as.mcmc(.x$samples)) %>% as.mcmc.list()
codaSamp <- codaSamp[,!grepl('state',colnames(codaSamp[[1]]))]

## Obtain WAIC value
waic <- chain_output %>% map_dbl(~.x$WAIC)
waic

## Save output
fit1 <- list(model = myCode, codaSamp = codaSamp, waic = waic, dur = dur)
# write_rds(fit1, 'out_fit1.rds', compress = 'xz')

## Check output
summary(codaSamp)

library(MCMCvis)
MCMCsummary(codaSamp, params = c('B.age','B.veg','B.dens','B.densVeg'), n.eff = TRUE, round = 3)
MCMCsummary(codaSamp, params = c('year.p','mean.p','sd.p'), n.eff = TRUE, round = 3)
MCMCsummary(codaSamp, params = c('Sigma.raw'), n.eff = TRUE, round = 3)
MCMCsummary(codaSamp, params = c('ageM'), n.eff = TRUE, round = 3)

# Assess MCMC convergence visually using traceplots
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

# Assess MCMC convergence formally using the Gelman-Rubin diagnostic (Gelman and Rubin 1992)
gelman.diag(codaSamp[,paste0('B.age[',1:nAge,']')])
gelman.diag(codaSamp[,paste0('B.veg[',1:nAge,']')])
gelman.diag(codaSamp[,paste0('B.dens[',1:nAge,']')])
gelman.diag(codaSamp[,paste0('B.densVeg[',1:nAge,']')])

gelman.diag(codaSamp[,paste0('year.p[',1:ntimes,']')])
gelman.diag(codaSamp[,'mean.p'])
gelman.diag(codaSamp[,'sd.p'])

gelman.diag(codaSamp[,paste0('Sigma.raw[',1:nAge,', ',1:nAge,']')])
gelman.diag(codaSamp[,paste0('ageM[',1:nb.noAge,']')])

# Check effective sample sizes
effectiveSize(codaSamp[,paste0('B.age[',1:nAge,']')])
effectiveSize(codaSamp[,paste0('B.veg[',1:nAge,']')])
effectiveSize(codaSamp[,paste0('B.dens[',1:nAge,']')])
effectiveSize(codaSamp[,paste0('B.densVeg[',1:nAge,']')])

effectiveSize(codaSamp[,paste0('year.p[',1:ntimes,']')])
effectiveSize(codaSamp[,'mean.p'])
effectiveSize(codaSamp[,'sd.p'])

effectiveSize(codaSamp[,paste0('Sigma.raw[',1:nAge,', ',1:nAge,']')])
effectiveSize(codaSamp[,paste0('ageM[',1:nb.noAge,']')])

## Visualize results
betaz = codaSamp %>% map(as.data.frame) %>% bind_rows()

# Check for correlations among fixed effects
tmp <- codaSamp[,grepl('B.',colnames(codaSamp[[1]]))] %>% map(as.data.frame) %>% bind_rows()
corrplot::corrplot(cor(tmp,use = 'p'))

# Check random effects among demographic rates
# Check variance-correlation matrix, with Sigma.raw on diagonal
VarCorrMatrix <- array(NA,dim = c(myconst$nAge,myconst$nAge,nrow(betaz)))

for (i in 1:myconst$nAge){
  VarCorrMatrix[i,i,] <- betaz[,paste0('xi[',i,']')]*
    sqrt(betaz[,paste0('Sigma.raw[',i,', ',i,']')])
} #i

for (j in 1:(myconst$nAge-1)){
  for (i in (j+1):myconst$nAge){
    VarCorrMatrix[j,i,] <- ( betaz[,paste0('Sigma.raw[',i,', ',j,']')])/
      sqrt(betaz[,paste0('Sigma.raw[',j,', ',j,']')]*
             betaz[,paste0('Sigma.raw[',i,', ',i,']')])    # Calculate correlation coefficients
  }
}

round(apply(VarCorrMatrix, 1:2, mean, na.rm = T), 2)                          # Extract means
round(apply(VarCorrMatrix, 1:2, quantile, prob = 0.025, na.rm = T), 2)        # and credible intervals
round(apply(VarCorrMatrix, 1:2, quantile, prob = 0.975, na.rm = T), 2)

# Calculate mean survival probabilities, by age class and year
df = expand.grid(age = 1:22, year = 1:12)
s.pred = matrix(NA, nrow = nrow(df), ncol = nrow(betaz))

for(i in 1:nrow(df)){
  s.pred[i,] <- betaz[,paste0('B.age[',mydata$SageC[df$age[i]],']')] +
    betaz[,paste0('gamma[',df$year[i],', ', mydata$SageC[df$age[i]],']')]
  df$ageC[i] <- mydata$SageC[df$age[i]]
}

df$s = inv.logit(apply(s.pred,1,mean))                     # Extract means
df$s.cil = inv.logit(apply(s.pred,1,quantile,0.025))       # and credible intervals
df$s.cih = inv.logit(apply(s.pred,1,quantile,0.975))

# Plot main results
df <- df %>% mutate(ageC = as.factor(ageC))

library(ggplot2)
library(scales)
plot <- ggplot(df, aes(x = year, y = s)) +
  geom_path(aes(colour = ageC)) +
  labs(x = "Year", y = "Survival", colour = "Age class") +
  scale_x_continuous(breaks = pretty_breaks()) +
  scale_y_continuous(breaks = pretty_breaks(), limits = c(0,1))
plot

