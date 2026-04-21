# 9 April 2025
# Run survival model

## Set up ----------------------------------------------------------------------

# set toggles
testRun <- TRUE
parallelRun <- TRUE
envEffectsS <- TRUE
ageClasses <- 6

# load packages
library(tidyverse)
library(lubridate)
library(beepr)
library(here)
library(boot)
library(coda)
library(nimble)
library(parallel)
library(nimbleEcology)

# load data
source('wrangleData_en.R')
enData <- wrangleData_en(dens.data = "data/abundanceData_Proteus.csv",
                         veg.data  = "data/biomass data April 2009 - Jan 2025_updated Feb2025.xlsx",
                         wea.data  = "data/Prom_Weather_2008-2023_updated Jan2025 RB.xlsx",
                         wind.data = "data/POWER_Point_Daily_20080101_20241231_10M.csv",
                         obs.data  = "data/PromObs_2008-2023.xlsx",
                         list      = "data/PromlistAllOct24.xlsx")

source('wrangleData_sv.R')
svData <- wrangleData_sv(surv.data = "data/PromSurvivalOct24.xlsx",
                         yafs.data = "data/RSmainRB_Mar25.xlsx",
                         ageClasses = ageClasses, known.age = TRUE)

# center & scale data
sc <- function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)

# center & scale density
dens  <- as.numeric(sc(enData$dens))
densE <- as.numeric(enData$densE/sd(enData$dens, na.rm = T))

# create Nimble lists
myData  <- list(obs    = svData$obs,
                # state  = svData$state,
                age.S  = svData$age.S,
                ageC.S = svData$ageC.S,
                
                dens  = dens,
                densE = densE,
                veg   = enData$veg,
                vegE  = enData$vegE,
                win   = enData$win)

myConst <- list(nID.S   = svData$nID.S,
                nYear   = svData$nYear,
                nAgeC.S = svData$nAgeC.S,
                
                first = svData$first,
                last  = svData$last,
                W     = diag(svData$nAgeC.S),
                DF    = svData$nAgeC.S,
                
                noDens  = enData$noDens,
                noVeg   = enData$noVeg,
                noWin   = enData$noWin,
                nNoDens = enData$nNoDens,
                nNoVeg  = enData$nNoVeg,
                nNoWin  = enData$nNoWin,
                
                envEffectsS = envEffectsS,
                ageClasses  = ageClasses
                # Si = rep(0, svData$nYear-1)
                )

# list2env(myData, envir = .GlobalEnv)
# list2env(myConst, envir = .GlobalEnv)

# # checks
# sapply(1:nrow(myData$state), function (i) myData$state[i, myConst$first[i]]) %>% table(useNA = 'a') # should all be 1
# sapply(1:nrow(myData$state), function (i) myData$state[i, myConst$last[i]]) %>% table(useNA = 'a')  # should be mostly 0s & NAs
# sapply(1:nrow(myData$state), function (i) myData$age[i, myConst$first[i]]) %>% table(useNA = 'a')   # should all be >= 1
# table(myConst$first >= myConst$last, useNA = 'a')                                                   # should all be F

# nimbleFunction to replace dCJS_vv
dCJS_fun <- nimbleFunction(
  run = function(x = double(2),
                 probSurvive = double(2),
                 probCapture = double(2),
                 len = integer(0),
                 log = logical(0)){
    
    returnType(double(0))
    
    # coerce matrix slices into true vectors
    vecX <- c(x)
    vecS <- c(probSurvive)
    vecC <- c(probCapture)
    
    if(len <= 0) return(0.0)
    
    # forward algorithm for CJS conditional on first capture at x[1]
    # compute joint prob of the observed history up to last occasion
    # alpha_alive: prob(alive & history up to current occasion)
    # alpha_dead: prob(dead & history up to current occasion)
    alpha_alive <- 1.0  # alive at first capture
    alpha_dead  <- 0.0
    
    for(j in 2:len){
      surv <- vecS[j-1] # survival prob between j-1 & j
      capt <- vecC[j]   # detection prob at occasion j
      
      
      if(vecX[j] == 1){
        alpha_alive_new <- alpha_alive * surv * capt            # if observed at j
      }else{
        alpha_alive_new <- alpha_alive * surv * (1.0 - capt)    # if not observed at j
      }
      alpha_dead_new <- alpha_dead + alpha_alive * (1.0 - surv) # if died between j-1 & j
      
      alpha_alive <- alpha_alive_new
      alpha_dead <- alpha_dead_new
    }
    
    prob <- alpha_alive + alpha_dead
    if(prob <= 0.0){
      return(-1e12)
    }
    return(log(prob))
  }
)

# Register this nimbleFunction as a distribution named 'dCJS_fun'
registerDistributions(list(
  dCJS_fun = list(BUGSdist = "dCJS_fun(x, probSurvive, probCapture, len)",
                  Rdist    = "dCJS_fun(x, probSurvive, probCapture, len, log)",
                  types    = c("x = double(2)",
                               "probSurvive = double(2)",
                               "probCapture = double(2)",
                               "len = integer(0)",
                               "log = logical(0)"))))


## Model -----------------------------------------------------------------------

myCode = nimbleCode({
  
  # #### Likelihood ####
  # for(i in 1:nID.S){
  #   for(t in (first[i] + 1):last[i]){
  #     # state process
  #     state[i, t] ~ dbern(Mu.Sp[i, t])
  #     Mu.Sp[i, t] <- S[ageC.S[age.S[i, t-1]], t-1] * state[i, t-1]
  #     
  #     # observation process
  #     obs[i, t] ~ dbern(Mu.Op[i, t])
  #     Mu.Op[i, t] <- O[t] * state[i, t]
  #   }
  # }
  
  # Using nimbleEcology
  for(i in 1:nID.S){
    
    # # ATTEMPT 1
    # Has everyone on the full timescale of 1:nYear
    # Does not work because most first[i] is often not at t = 1
    # Nimble therefore does not accept this because obs is often not 1 at t = 1
    
    # for(t in 1:(nYear-1)){
    #   Si[t] <- S[ageC.S[age.S[i, t]], t] * (t >= first[i] & t <= (last[i]-1))
    # }
    # 
    # obs[i, 1:nYear] ~ dCJS_vv(probSurvive = Si[1:(nYear-1)],
    #                           probCapture = O[1:nYear],
    #                           len = nYear)
    
    # # ATTEMPT 2
    # Has every roo on its own timescale of first:last
    # Does not work because probSurvive requires a true 1-D vector
    # & simply refuses to see Si[first[i]:(last[i]-1)] as a vector
    
    # for(t in first[i]:(last[i]-1)){
    #   Si[t] <- S[ageC.S[age.S[i, t]], t]
    # }
    # 
    # obs[i, first[i]:last[i]] ~ dCJS_vv(probSurvive = Si[first[i]:(last[i]-1)],
    #                                    probCapture = O[first[i]:last[i]],
    #                                    len = last[i] - first[i] + 1)
    
    # ATTEMPT 3
    # Customized nimbleFunction: a Hail Mary, entirely chatGPT's concoction
    # Should be replicating dCJS_vv, & coercing matrix slices to vectors
    # Does not work: "bad parameters for distribution dCJS_fun" error
    
    for(t in first[i]:(last[i]-1)){
      Si[i, t] <- S[ageC.S[age.S[i, t]], t]
    }
    
    obs[i, first[i]:last[i]] ~ dCJS_fun(probSurvive = Si[i, first[i]:(last[i]-1)], # 1D matrix slice
                                        probCapture = O[first[i]:last[i]],
                                        len = last[i] - first[i] + 1)
  }
  
  #### Constraints ####
  # survival function
  for(a in 1:nAgeC.S){                               
    for(t in 1:(nYear-1)){
      if(envEffectsS){
        logit(S[a, t]) <- logit(Mu.S[a]) +
          BetaD.S * dens[t] +
          BetaV.S * veg[t] +
          BetaW.S * win[t] +
          # Gamma.S[t, a]
          EpsilonT.S[t]
      }else{
        logit(S[a, t]) <- logit(Mu.S[a]) +
          # Gamma.S[t, a]
          EpsilonT.S[t]
      }
    }
  }
  
  # observation function
  for(t in 1:nYear){
    EpsilonT.O[t] ~ dnorm(0, sd = SigmaT.O)
    logit(O[t]) <- logit(Mu.O) + EpsilonT.O[t]
  }
  
  # missing environment
  if(envEffectsS){
    for(m in 1:nNoDens) dens[noDens[m]] ~ dnorm(0, sd = 1)
    for(m in 1:nNoVeg) veg[noVeg[m]] ~ dnorm(0, sd = 1)
    for(m in 1:nNoWin) win[noWin[m]] ~ dnorm(0, sd = 1)
  }
  
  #### Priors ####
  # for intercepts
  for(a in 1:nAgeC.S){
    Mu.S[a] ~ dunif(0.01, 0.99)
  }
  
  # # for age-dependent fixed effects
  # if(envEffectsS){
  #   for(a in 1:nAgeC.S){
  #     BetaD.S[a] ~ dunif(-5, 5)
  #     BetaV.S[a] ~ dunif(-5, 5)
  #     BetaW.S[a] ~ dunif(-5, 5)
  #   }
  # }
  
  # for age-independent fixed effects
  if(envEffectsS){
    BetaD.S ~ dunif(-5, 5)
    BetaV.S ~ dunif(-5, 5)
    BetaW.S ~ dunif(-5, 5)
  }
  
  # # for age-dependent random effects
  # # variance-covariance matrix
  # for(a in 1:nAgeC.S){
  #   zero[a] <- 0
  #   Xi.S[a] ~ dunif(0, 2)
  # }
  # 
  # for(t in 1:(nYear-1)){
  #   Epsilon.S[t, 1:nAgeC.S] ~ dmnorm(zero[1:nAgeC.S], Tau.S[1:nAgeC.S, 1:nAgeC.S])
  #   for(a in 1:nAgeC.S){
  #     Gamma.S[t, a] <- Xi.S[a] * Epsilon.S[t, a]
  #   }
  # }
  # 
  # # precision matrix
  # Tau.S[1:nAgeC.S, 1:nAgeC.S] ~ dwish(W[1:nAgeC.S, 1:nAgeC.S], DF)
  # Sigma.S[1:nAgeC.S, 1:nAgeC.S] <- inverse(Tau.S[1:nAgeC.S, 1:nAgeC.S])
  
  # for age-independent random effects
  for(t in 1:(nYear-1)){
    XiT.S[t] ~ dnorm(0, sd = 1) # latent standard normal
    EpsilonT.S[t] <- SigmaT.S * XiT.S[t] # actual random effect
  }
  SigmaT.S ~ dunif(0, 10) # scale of the random effect
  
  # for observation
  Mu.O ~ dunif(0.01, 0.99) # or dunif(0, 1)
  SigmaT.O ~ dunif(0.01, 10) # or dunif(0, 10)
  
})


## Assemble --------------------------------------------------------------------

nchains   <- 4
seedMod   <- 1:nchains
seedInits <- 1

# Tau.S = diag(myConst$nAgeC.S) + rnorm(myConst$nAgeC.S^2, 0, 0.1)
# Tau.S = (Tau.S + t(Tau.S))/2

# # alternative init Tau.S that guarantees positive-definite values
# # (apparently potentially problematic the way I had it above)
# A <- matrix(rnorm(nAgeC.S^2, 0, 0.1), nAgeC.S, nAgeC.S)
# Tau.S <- crossprod(A) + diag(nAgeC.S)  # positive-definite

XiT.S = rnorm(myConst$nYear-1, 0, 1)
SigmaT.S = runif(1, .5, 2)
EpsilonT.S = XiT.S * SigmaT.S

set.seed(seedInits)
myInits <- list()
for(c in 1:nchains){
myInits[[c]] <- list(
  dens = round(ifelse(is.na(myData$dens), rnorm(myConst$nYear, 0, .1), myData$dens), 4),
  veg  = round(ifelse(is.na(myData$veg), rnorm(myConst$nYear, 0, .1), myData$veg), 4),
  win  = round(ifelse(is.na(myData$win), rnorm(myConst$nYear, 0, .1), myData$win), 4),
  
  Mu.S = rep(runif(myConst$nAgeC.S, 0.1, 0.9)),
  
  O = runif(myConst$nYear, 0.1, 0.9),
  Mu.O = runif(1, 0.1, 0.9),
  EpsilonT.O = rnorm(myConst$nYear, 0, 0.2),
  SigmaT.O = runif(1, 0.01, 2),
  
  # # for age-dependent random effect
  # Xi.S = rnorm(myConst$nAgeC.S, 1, 0.1),
  # Epsilon.S = matrix(rnorm((myConst$nYear-1) * myConst$nAgeC.S, 0, 0.1),
  #                    nrow = myConst$nYear-1, ncol = myConst$nAgeC.S),
  # Tau.S = Tau.S
  
  # for age-independent random effect
  XiT.S = XiT.S,
  SigmaT.S = SigmaT.S,
  EpsilonT.S = EpsilonT.S
)}

# select parameters to monitors
params <- c('S', 'Mu.S',
            # 'Gamma.S', 'Sigma.S',            # random effects (correlated)
            'EpsilonT.S', 'SigmaT.S',          # random effects (uncorrelated)
            'Mu.O', 'EpsilonT.O', 'SigmaT.O'
            # 'dens.true', 'veg.true', 'win.true'
            )

if(envEffectsS){params <- c(params, 'BetaD.S', 'BetaV.S', 'BetaW.S')}

# select MCMC settings
if(testRun){
  nthin   <- 1
  nburnin <- 0
  niter   <- 10
}else{
  nthin   <- 8
  nburnin <- 40000
  niter   <- nburnin + 1000*nthin
}


## Run model -------------------------------------------------------------------

if(parallelRun){
  # function to run one chain inside cluster
  runChain <- function(chainID, code, data, const, inits, params,
                       niter, nburnin, nthin, seed){
    
    library(nimble)
    library(nimbleEcology)
    
    registerDistributions(list(
      dCJS_fun = list(BUGSdist = "dCJS_fun(x, probSurvive, probCapture, len)",
                      Rdist    = "dCJS_fun(x, probSurvive, probCapture, len, log)",
                      types    = c("x = double(2)",
                                   "probSurvive = double(2)",
                                   "probCapture = double(2)",
                                   "len = integer(0)",
                                   "log = logical(0)"))))
    
    set.seed(seed)
    inits <- myInits[[chainID]]
    
    model <- nimbleModel(code = code, data = data, constants = const, inits = inits)
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
                                "seedMod", "runChain",
                                "dCJS_fun"))
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
saveRDS(out.mcmc, 'results/CJS_dCJSvv.rds', compress = 'xz')


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
MCMCsummary(out.mcmc, params = c('BetaD.S', 'BetaV.S', 'BetaW.S'), n.eff = TRUE, round = 2)
MCMCsummary(out.mcmc, params = c('Mu.O', 'EpsilonT.O', 'SigmaT.O'), n.eff = TRUE, round = 2)
MCMCsummary(out.mcmc, params = c('Sigma.S'), n.eff = TRUE, round = 2)

# chainplots
MCMCtrace(out.mcmc, params = c('S', 'BetaA.S'), pdf = FALSE)
MCMCtrace(out.mcmc, params = c('BetaD.S', 'BetaV.S', 'BetaW.S'), pdf = FALSE)
MCMCtrace(out.mcmc, params = c('Mu.O', 'EpsilonT.O', 'SigmaT.O'), pdf = FALSE)
MCMCtrace(out.mcmc, params = c('Sigma.S'), pdf = FALSE)


## Plots -----------------------------------------------------------------------

out.dat <- out.mcmc %>% map(as.data.frame) %>% bind_rows()

# check for correlations among fixed effects
par(mfrow = c(1,1))
corrplot(cor(out.mcmc[, grepl('Beta', colnames(out.mcmc[[1]]))] %>%
               map(as.data.frame) %>% bind_rows(), use = 'p'))

# check random effects among demographic rates
# check variance-correlation matrix, with Sigma.S on diagonal
varCorrMatrix <- array(NA, dim = c(myConst$nAgeC.S, myConst$nAgeC.S, nrow(out.dat)))

for(i in 1:myConst$nAgeC.S){
  varCorrMatrix[i,i,] <- out.dat[, paste0('Xi.S[', i,']')]*
    sqrt(out.dat[, paste0('Sigma.S[', i,', ', i,']')])
}

for(j in 1:(myConst$nAgeC.S-1)){
  for(i in (j+1):myConst$nAgeC.S){
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
  S.pred[i, ] <- out.dat[, paste0('BetaA.S[', myData$ageC.S[df$age[i]], ']')] +
    out.dat[, paste0('Gamma.S[', df$year[i], ', ', myData$ageC.S[df$age[i]], ']')]
  df$ageC.S[i] <- myData$ageC.S[df$age[i]]
}

df$S = inv.logit(apply(S.pred, 1, mean))
df$sLCI = inv.logit(apply(S.pred, 1, quantile, 0.025))
df$sUCI = inv.logit(apply(S.pred, 1, quantile, 0.975))

# plot main results
df %>%
  mutate(ageC.S = as.factor(ageC.S)) %>% 
  ggplot(aes(x = year, y = S)) +
  geom_ribbon(aes(ymin = sLCI, ymax = sUCI, fill = ageC.S), alpha = 0.1) +
  geom_line(aes(colour = ageC.S), linewidth = 1, show.legend = F) +
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

