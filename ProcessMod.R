# 8 April 2025
# Process model for roo IPM
# 1 census per year, on Sept 1 ish

## Set up ----------------------------------------------------------------------

# set toggles
testRun <- FALSE

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

source("wrangleData_rs.R")
rsData <- wrangleData_rs(rs.data = "data/RSmainRB_Mar25.xlsx",
                         obs.data = "data/PromObs_2008-2019.xlsx",
                         known.age = TRUE, cum.surv = TRUE, surv.sep1 = TRUE)

# create Nimble lists
myData  <- list(obs = svData$obs,
                state = svData$state,
                age = svData$age,
                ageC = svData$ageC,
                dens = enData$dens,
                densE = enData$densE,
                veg = enData$veg,
                vegE = enData$vegE)

myConst <- list(nYear = rsData$nYear,
                nIDs = svData$nID,
                nIDr = rsData$nID,
                nAge = rsData$nAge,
                nAgeC = rsData$nAgeC,
                noAge = svData$noAge,
                nNoAge = svData$nNoAge,
                nNoVeg = enData$nNoVeg,
                first = svData$first,
                last = svData$last,
                W = svData$W,
                DF = svData$DF)


## Assemble --------------------------------------------------------------------

source("writeCode.R")
myCode <- writeCode()

nchains   <- 3
seedMod   <- 1:nchains
seedInits <- 1

# assign initial values
source("simulateInits.R")
set.seed(seedInits)
myInits <- list()
for(c in 1:nchains){
  myInits[[c]] <- simulateInits(
    # n = rsData$n,
    nIDs = myConst$nIDs,
    # nIDr = myConst$nIDr,
    nYear = myConst$nYear,
    nAge = myConst$nAge,
    nAgeC = myConst$nAgeC,
    
    age = myData$age,
    dens = myData$dens,
    veg = myData$veg,
    # win = myData$win,
    
    nNoAge = myConst$nNoAge,
    # nNoDens = myConst$nNoDens,
    nNoVeg = myConst$nNoVeg,
    # nNoWin = myConst$nNoWin
    )
}

# select parameters to monitors
params = c(
  # Survival model
  'dens.hat', 'veg.hat', # 'ageM',            # latent states
  'BetaA.sv', 'BetaD.sv', 'BetaV.sv',         # covariate effects
  'Mu.ob', 'Epsilon.ob', 'Sigma.ob',          # observation parameters
  'Gamma.sv', 'Xi.sv', 'Sigma.sv',            # random effects
  
  # Population model
  'sv', 'b', 'svPY', 'svYAF', 'svSA', 'svAD', # yearly vital rates
  'nYAF', 'nSA', 'nAD', 'nTOT')               # population sizes

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

# MCMC output
out.mcmc <- as.mcmc.list(samples)

# save output
fit <- list(model = myCode, out.mcmc = out.mcmc, dur = dur)
# write_rds(fit, 'results/IPM_CJS.rds', compress = 'xz')


## Results ---------------------------------------------------------------------

library(coda)
library(MCMCvis)
library(corrplot)
library(ggplot2)
library(scales)

# summary(out.mcmc) # cannot handle NAs

# # find parameters generating NAs
# for(i in 1:ncol(out.mcmc[[1]])){
#   if(any(is.na(out.mcmc[[1]][,i]))){
#     message(paste0(colnames(out.mcmc[[1]])[i]))
#     print(out.mcmc[[1]][1:3,i])
#   }
# }

# summaries
MCMCsummary(out.mcmc, params = c('BetaA.sv', 'BetaD.sv', 'BetaV.sv'), n.eff = TRUE, round = 2)
MCMCsummary(out.mcmc, params = c('Mu.ob', 'Epsilon.ob', 'Sigma.ob'), n.eff = TRUE, round = 2)
MCMCsummary(out.mcmc, params = c('Sigma.sv'), n.eff = TRUE, round = 2)

MCMCsummary(out.mcmc, params = c('sv', 'b', 'svPY'), n.eff = TRUE, round = 2)
MCMCsummary(out.mcmc, params = c('svYAF', 'svSA', 'svAD'), n.eff = TRUE, round = 2)
MCMCsummary(out.mcmc, params = c('nYAF', 'nSA', 'nAD', 'nTOT'), n.eff = TRUE, round = 2)

# chainplots
MCMCtrace(out.mcmc, params = c('BetaA.sv', 'BetaD.sv', 'BetaV.sv'), pdf = FALSE)
MCMCtrace(out.mcmc, params = c('Epsilon.ob', 'Sigma.ob'), pdf = FALSE)
MCMCtrace(out.mcmc, params = c('Sigma.sv'), pdf = FALSE)

MCMCtrace(out.mcmc, params = c('sv', 'b', 'svPY'), pdf = FALSE)
MCMCtrace(out.mcmc, params = c('svYAF', 'svSA', 'svAD'), pdf = FALSE)
MCMCtrace(out.mcmc, params = c('nYAF', 'nSA', 'nAD', 'nTOT'), pdf = FALSE)


## Plot population model -------------------------------------------------------

nYear <- myConst$nYear
nAge  <- myConst$nAge
nAgeC <- myConst$nAgeC

# posterior samples
out.mat <- as.matrix(samples)

# parameters to include
table.params <- c(
  paste0('nYAF[', 1:nYear, ']'),
  paste0('nSA[', rep(1:2, each = nYear), ', ', rep(1:nYear, times = 2), ']'),
  paste0('nAD[', rep(1:(nAge+2), each = nYear), ', ', rep(1:nYear, times = (nAge+2)), ']'))

# table.params <- list(
#   nYAF = c(paste0('nYAF[', 1:nYear, ']')),
#   nSA  = c(paste0('nSA[', rep(1:2, each = nYear), ', ', rep(1:nYear, times = 2), ']')),
#   nAD  = c(paste0('nAD[', rep(1:nAge, each = nYear), ', ', rep(1:nYear, times = nAge), ']')))

# table of posterior summaries
post.table <- data.frame(Parameter = table.params, Estimate = NA)

for(i in 1:length(table.params)){
  est <- out.mat[, table.params[i]]
  post.table$Estimate[i] <- paste0(round(median(est, na.rm = T), digits = 2), ' [',
                                   round(quantile(est, 0.025, na.rm = T), digits = 2), ', ',
                                   round(quantile(est, 0.975, na.rm = T), digits = 2), ']')
}

# plot results
# nTOT <- grep("^nTOT\\[", colnames(samples)); nTOT
# nYAF <- grep("^nYAF\\[", colnames(samples)); nYAF
# nSA <- grep("^nSA\\[", colnames(samples)); nSA
# nAD <- grep("^nAD\\[", colnames(samples)); nAD

nTOT <- grep("^nTOT\\[", colnames(out.mcmc[[1]])); nTOT
nYAF <- grep("^nYAF\\[", colnames(out.mcmc[[1]])); nYAF
nSA <- grep("^nSA\\[", colnames(out.mcmc[[1]])); nSA
nAD <- grep("^nAD\\[", colnames(out.mcmc[[1]])); nAD

var <- nTOT

df <- data.frame(
  Year = 1:length(var),
  Mean = apply(out.mat[, var, drop = FALSE], 2, mean, na.rm = TRUE),
  Lower = apply(out.mat[, var, drop = FALSE], 2, quantile, probs = 0.025, na.rm = TRUE),
  Upper = apply(out.mat[, var, drop = FALSE], 2, quantile, probs = 0.975, na.rm = TRUE)
)

# population plot
ggplot(df, aes(x = Year, y = Mean)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "#C398B7", alpha = 0.4) +
  geom_line(color = "#673C5B", linewidth = 1) +
  ylab("Parameter value") +
  xlab("Year") +
  theme_bw()

# ggsave("figures/IPM.jpeg", scale = 1, width = 18.0, height = 9.0, units = c("cm"), dpi = 600)


## Plots survival model --------------------------------------------------------

out.dat <- out.mcmc %>% map(as.data.frame) %>% bind_rows()

# check for correlations among fixed effects
par(mfrow = c(1,1))
corrplot(cor(out.mcmc[, grepl('B.', colnames(out.mcmc[[1]]))] %>%
               map(as.data.frame) %>% bind_rows(), use = 'p'))

# check random effects among demographic rates
# check variance-correlation matrix, with Sigma.sv on diagonal
varCorrMatrix <- array(NA, dim = c(myConst$nAgeC, myConst$nAgeC, nrow(out.dat)))

for(i in 1:myConst$nAgeC){
  varCorrMatrix[i,i,] <- out.dat[, paste0('Xi.sv[', i,']')]*
    sqrt(out.dat[, paste0('Sigma.sv[', i,', ', i,']')])
}

for(j in 1:(myConst$nAgeC-1)){
  for(i in (j+1):myConst$nAgeC){
    varCorrMatrix[j, i, ] <- (out.dat[, paste0('Sigma.sv[', i, ', ', j, ']')])/
      sqrt(out.dat[, paste0('Sigma.sv[', j, ', ', j, ']')]*
             out.dat[, paste0('Sigma.sv[', i, ', ', i, ']')])
  }
}

round(apply(varCorrMatrix, 1:2, mean, na.rm = T), 2)
round(apply(varCorrMatrix, 1:2, quantile, prob = 0.025, na.rm = T), 2)
round(apply(varCorrMatrix, 1:2, quantile, prob = 0.975, na.rm = T), 2)

# calculate survival probabilities
df <- expand.grid(age = 1:22, year = 1:16)
sv.pred <- matrix(NA, nrow = nrow(df), ncol = nrow(out.dat))

for(i in 1:nrow(df)){
  sv.pred[i, ] <- out.dat[, paste0('BetaA.sv[', myData$ageC[df$age[i]], ']')] +
    out.dat[, paste0('Gamma.sv[', df$year[i], ', ', myData$ageC[df$age[i]], ']')]
  df$ageC[i] <- myData$ageC[df$age[i]]
}

df$sv = inv.logit(apply(sv.pred, 1, mean))
df$svLCI = inv.logit(apply(sv.pred, 1, quantile, 0.025))
df$svUCI = inv.logit(apply(sv.pred, 1, quantile, 0.975))

# plot main results
df %>%
  mutate(ageC = as.factor(ageC)) %>% 
  ggplot(aes(x = year, y = sv)) +
  geom_ribbon(aes(ymin = svLCI, ymax = svUCI, fill = ageC), alpha = 0.2) +
  geom_line(aes(colour = ageC), linewidth = 1, show.legend = F) +
  labs(x = "Year", y = "Survival", fill = "Age class") +
  # scale_x_continuous(breaks = pretty_breaks()) +
  scale_y_continuous(breaks = pretty_breaks()) +
  theme_bw()

# ggsave("figures/IPM_CJS.jpeg", scale = 1, width = 18.0, height = 9.0, units = c("cm"), dpi = 600)

