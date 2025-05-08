# 22 April 2025
# Run reproductive success model

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
source("wrangleData_en.R")
enData <- wrangleData_en(dens.data = "data/abundanceData_Proteus.csv",
                         veg.data  = "data/biomass data April 2009 - Jan 2025_updated Feb2025.xlsx",
                         wea.data  = "data/Prom_Weather_2008-2023_updated Jan2025 RB.xlsx",
                         wind.data = "data/POWER_Point_Daily_20080101_20241231_10M.csv",
                         obs.data  = "data/PromObs_2008-2019.xlsx",
                         list      = "data/PromlistAllOct24.xlsx")

source("wrangleData_rs.R")
rsData <- wrangleData_rs(rs.data = "data/RSmainRB_Mar25.xlsx",
                         obs.data = "data/PromObs_2008-2019.xlsx",
                         known.age = TRUE, cum.surv = FALSE)

# to play around with age classes!
# ageC <- c(1,2,2,3,3,3,3,4,4,4, rep(5,30)) # default
ageC <- c(seq(from = 1, to = 20, by = 1), rep(20, times = 20)); ageC

# create Nimble lists
myData <-  list(B       = rsData$B,
                R       = rsData$survS1,
                id.R    = rsData$id.R,
                year.R  = rsData$year.R,
                age.R   = rsData$age.R,
                dens    = enData$dens,
                densE   = enData$densE,
                veg     = enData$veg,
                vegE    = enData$vegE,
                win     = enData$win)

myConst <- list(nR      = rsData$nR,
                nID.R   = rsData$nID.R,
                nYear   = rsData$nYear,
                nAge    = rsData$nAge,
                nAgeC   = max(ageC),
                nNoDens = enData$nNoDens,
                nNoVeg  = enData$nNoVeg,
                nNoWin  = enData$nNoWin)


## Model -----------------------------------------------------------------------

myCode = nimbleCode({
  
  #### Likelihood & constraints ####
  # yearly birth rate
  for(x in 1:nR){
    B[x] ~ dbern(Bi[x])
    logit(Bi[x]) <- logit(Mu.B) +
      EpsilonT.B[year.R[x]]
  }
  
  for(t in 1:(nYear-1)){
    logit(Bt[t]) <- logit(Mu.B) + EpsilonT.B[t]
  }
  
  # individual RS function
  Mu.Ri[1] <- 0
  for(x in 1:nR){
    R[x] ~ dbern(Ri[x])
    logit(Ri[x]) <- logit(Mu.Ri[age.R[x]]) +
      BetaD.R * dens.hat[year.R[x]] +
      BetaV.R * veg.hat[year.R[x]] +
      BetaW.R * win.hat[year.R[x]] +
      EpsilonI.Ri[id.R[x]] +
      EpsilonT.Ri[year.R[x]]
  }
  
  # age-specific RS function
  # use parameters estimated from individual data above
  # to predict age-specific reproductive success (Ra) here!
  Mu.Ra[1] <- 0
  for(a in 1:nAge){
    for(t in 1:(nYear-1)){
      logit(Ra[a, t]) <- logit(Mu.Ra[a]) +
        BetaD.R * dens.hat[t] +
        BetaV.R * veg.hat[t] +
        BetaW.R * win.hat[t] +
        EpsilonT.Ra[t]
    }
  }
  
  # missing environment
  for(t in 1:(nYear-1)){
    dens.hat[t] ~ dnorm(dens[t], sd = densE[t])
    veg.hat[t]  ~ dnorm(veg[t], sd = vegE[t])
    win.hat[t]  ~ dnorm(win[t], sd = 1)
  }
  
  for(m in 1:nNoVeg){
    veg[m] ~ dnorm(0, sd = 2)
  }
  
  for(m in 1:nNoDens){
    dens[m] ~ dnorm(0, sd = 2)
  }
  
  for(m in 1:nNoWin){
    win[m] ~ dnorm(0, sd = 2)
  }
  
  ##### Priors ####
  # priors for fixed effects
  for(a in 2:nAge){
    Mu.Ri[a] ~ dunif(0, 1)
    Mu.Ra[a] ~ dunif(0, 1)
  }
  Mu.B ~ dunif(0, 1)
  
  BetaD.R ~ dunif(-2, 2) # could be dunif(-5, 5) if need be
  BetaV.R ~ dunif(-2, 2) # could be dunif(-5, 5) if need be
  BetaW.R ~ dunif(-2, 2) # could be dunif(-5, 5) if need be
  
  # priors for random effects
  for(i in 1:nID.R){
    EpsilonI.Ri[i] ~ dnorm(0, sd = SigmaI.Ri)
  }
  
  for(t in 1:(nYear-1)){
    EpsilonT.Ri[t] ~ dnorm(0, sd = SigmaT.Ri)
    EpsilonT.Ra[t] ~ dnorm(0, sd = SigmaT.Ra)
    EpsilonT.B[t] ~ dnorm(0, sd = SigmaT.B)
  }
  
  # priors for sigma
  SigmaI.Ri ~ dunif(0, 100)
  SigmaT.Ri ~ dunif(0, 100)
  SigmaT.Ra ~ dunif(0, 100)
  SigmaT.B ~ dunif(0, 100)
  
})


## Assemble --------------------------------------------------------------------

nchains   <- 4
seedMod   <- 1:nchains
seedInits <- 1

# assign initial values
source("simulateInits.R")
set.seed(seedInits)
myInits <- list()
for(c in 1:nchains){
  myInits[[c]] <- simulateInits(
    nR = rsData$nR,
    nID.R = myConst$nID.R,
    nYear = myConst$nYear,
    nAge = myConst$nAge,
    nAgeC = myConst$nAgeC,
    id.R = myData$id.R,
    year.R = myData$year.R,
    age.R = myData$age.R,
    dens = enData$dens,
    veg = enData$veg,
    win = enData$win,
    propF = enData$propF,,
    envEffectsR = TRUE,
    envEffectsS = TRUE
  )
}

# select parameters to monitor
params = c("Bt", "Ra",
           "Mu.B", "Mu.Ri", "Mu.Ra",
           "BetaD.R", "BetaV.R", "BetaW.R",
           "EpsilonT.B", "EpsilonI.Ri", "EpsilonT.Ri", "EpsilonT.Ra", 
           "SigmaT.B", "SigmaI.Ri", "SigmaT.Ri", "SigmaT.Ra",
           "dens.hat", "veg.hat", "win.hat"
)

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
# saveRDS(out.mcmc, 'results/RS_Age20.rds', compress = 'xz')


## Checks ----------------------------------------------------------------------

library(coda)
library(MCMCvis)
library(corrplot)
library(ggplot2)
library(scales)
library(patchwork)

# # load results
# out.mcmc <- readRDS('results/RS_Age20.rds')
# summary(out.mcmc) # cannot handle NAs

# summaries
MCMCsummary(out.mcmc, params = c('Bt', 'Ra'), n.eff = TRUE, round = 2)
MCMCsummary(out.mcmc, params = c('Mu.B', 'Mu.Ri', 'Mu.Ra'), n.eff = TRUE, round = 2)
MCMCsummary(out.mcmc, params = c('BetaD.R', 'BetaV.R', 'BetaW.R'), n.eff = TRUE, round = 2)
MCMCsummary(out.mcmc, params = c('SigmaT.B', 'SigmaI.Ri', 'SigmaT.Ri', 'SigmaT.Ra'), n.eff = TRUE, round = 2)

# chainplots
MCMCtrace(out.mcmc, params = c('Bt', 'Ra'), pdf = FALSE)
MCMCtrace(out.mcmc, params = c('Mu.B', 'Mu.Ri', 'Mu.Ra'), pdf = FALSE)
MCMCtrace(out.mcmc, params = c('BetaD.R', 'BetaV.R', 'BetaW.R'), pdf = FALSE)
MCMCtrace(out.mcmc, params = c('SigmaT.B', 'SigmaI.Ri', 'SigmaT.Ri', 'SigmaT.Ra'), pdf = FALSE)


## Plots -----------------------------------------------------------------------

out.dat <- out.mcmc %>% map(as.data.frame) %>% bind_rows()

# calculate reproductive probabilities
df <- expand.grid(age = 1:19, year = 1:16)
Bt.pred <- matrix(NA, nrow = nrow(df), ncol = nrow(out.dat))
Ra.pred <- matrix(NA, nrow = nrow(df), ncol = nrow(out.dat))

for(i in 1:nrow(df)){
  Bt.pred[i, ] <- qlogis(out.dat[, "Mu.B"]) +
    out.dat[, paste0("EpsilonT.B[", df$year[i], "]")]
}

df$Bt <- inv.logit(apply(Bt.pred, 1, mean))
df$bLCI <- inv.logit(apply(Bt.pred, 1, quantile, 0.025))
df$bUCI <- inv.logit(apply(Bt.pred, 1, quantile, 0.975))

for(i in 1:nrow(df)){
  Ra.pred[i, ] <- qlogis(out.dat[, paste0("Mu.Ra[", df$age[i], "]")]) +
    out.dat[, "BetaD.R"] * out.dat[, paste0("dens.hat[", df$year[i], "]")] +
    out.dat[, "BetaV.R"] * out.dat[, paste0("veg.hat[", df$year[i], "]")] +
    out.dat[, "BetaW.R"] * out.dat[, paste0("win.hat[", df$year[i], "]")] +
    out.dat[, paste0("EpsilonT.Ra[", df$year[i], "]")]
}

df$Ra <- inv.logit(apply(Ra.pred, 1, mean))
df$rLCI <- inv.logit(apply(Ra.pred, 1, quantile, 0.025))
df$rUCI <- inv.logit(apply(Ra.pred, 1, quantile, 0.975))

# plot main results
young <- df %>% 
  mutate(age = as.factor(age)) %>% 
  filter(age %in% c(2, 4, 6, 8, 10)) %>% 
  ggplot(aes(x = year, y = Ra)) +
  geom_ribbon(aes(ymin = rLCI, ymax = rUCI, fill = age), alpha = 0.1, show.legend = F) +
  geom_line(aes(colour = age), linewidth = 1) +
  labs(x = "Year", y = "Probability of a bean making it to YAF", colour = "Mom's age") +
  scale_y_continuous(breaks = pretty_breaks()) +
  theme_bw(); young
  
old <- df %>%
  mutate(age = as.factor(age-1)) %>% 
  filter(age %in% c(12, 14, 16, 18)) %>%
  ggplot(aes(x = year, y = Ra)) +
  geom_ribbon(aes(ymin = rLCI, ymax = rUCI, fill = age), alpha = 0.1, show.legend = F) +
  geom_line(aes(colour = age), linewidth = 1) +
  labs(x = "Year", y = "Probability of a bean making it to YAF", colour = "Mom's age") +
  scale_y_continuous(breaks = pretty_breaks()) +
  theme_bw(); old
  
young / old

# ggsave("figures/RS_AgeYr.jpeg", scale = 1, width = 18.0, height = 18.0, units = c("cm"), dpi = 600)

df %>% 
  group_by(year) %>% 
  mutate(Bt = mean(Bt),
         bUCI = mean(bUCI),
         bLCI = mean(bLCI)) %>% 
  ggplot(aes(x = year, y = Bt)) +
  geom_ribbon(aes(ymin = bLCI, ymax = bUCI), alpha = 0.1) +
  geom_line(linewidth = 1, show.legend = F) +
  scale_y_continuous(breaks = pretty_breaks()) +
  labs(x = "Year", y = "Birth rate") +
  ylim(0, 1) +
  theme_bw()

# ggsave("figures/RS_Bt.jpeg", scale = 1, width = 18.0, height = 9.0, units = c("cm"), dpi = 600)

df %>% 
  group_by(age) %>% 
  mutate(Ra = mean(Ra),
         rUCI = mean(rUCI),
         rLCI = mean(rLCI)) %>% 
  ggplot(aes(x = age, y = Ra)) +
  geom_ribbon(aes(ymin = rLCI, ymax = rUCI), alpha = 0.1) +
  geom_line(linewidth = 1, show.legend = F) +
  scale_y_continuous(breaks = pretty_breaks()) +
  labs(x = "Mom's age", y = "Probability of a bean making it to YAF") +
  theme_bw()

# ggsave("figures/RS_Ra.jpeg", scale = 1, width = 18.0, height = 9.0, units = c("cm"), dpi = 600)

