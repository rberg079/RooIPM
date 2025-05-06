# 8 April 2025
# Process model for roo IPM
# 1 census per year, on Sept 1 ish

## Set up ----------------------------------------------------------------------

# set toggles
testRun <- F

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

source('wrangleData_rs.R')
rsData <- wrangleData_rs(rs.data = "data/RSmainRB_Mar25.xlsx",
                         obs.data = "data/PromObs_2008-2019.xlsx",
                         known.age = TRUE, cum.surv = FALSE)

# create Nimble lists
myData  <- list(obs = svData$obs,
                state = svData$state,
                age.S = svData$age.S,
                ageC = svData$ageC,
                
                B = rsData$B,
                R = rsData$survS1,
                id.R = rsData$id.R,
                year.R = rsData$year.R,
                age.R = rsData$age.R,
                
                ab = enData$ab,
                abE = enData$abE,
                propF = enData$propF,
                dens = enData$dens,
                densE = enData$densE,
                veg = enData$veg,
                vegE = enData$vegE,
                win = enData$win)

myConst <- list(nR = rsData$nR,
                nID.S = svData$nID,
                nID.R = rsData$nID,
                nYear = svData$nYear,
                nAge = rsData$nAge,
                nAgeC = rsData$nAgeC,
                # noAge = svData$noAge,
                # nNoAge = svData$nNoAge,
                nNoDens = enData$nNoDens,
                nNoVeg = enData$nNoVeg,
                nNoWin = enData$nNoWin,
                nNoProp = enData$nNoProp,
                first = svData$first,
                last = svData$last,
                W = svData$W,
                DF = svData$DF)


## Assemble --------------------------------------------------------------------

source('writeCode.R')
myCode <- writeCode()

nchains   <- 4
seedMod   <- 1:nchains
seedInits <- 1

# assign initial values
source('simulateInits.R')
set.seed(seedInits)
myInits <- list()
for(c in 1:nchains){
  myInits[[c]] <- simulateInits(
    nR = rsData$nR,
    nID.S = myConst$nID.S,
    nID.R = myConst$nID.R,
    nYear = myConst$nYear,
    nAge = myConst$nAge,
    nAgeC = myConst$nAgeC,
    id.R = myData$id.R,
    year.R = myData$year.R,
    age.R = myData$age.R,
    dens = myData$dens,
    veg = myData$veg,
    win = myData$win,
    propF = myData$propF
    # nNoAge = myConst$nNoAge
    )
}

# select parameters to monitors
params = c(
  # Population model
  'S', 'Bt', 'Ra', 'sYAF', 'sSA', 'sAD',                      # yearly vital rates
  'nYAF', 'nSA', 'nAD', 'nTOT',                               # population sizes
  
  # Survival model
  'dens.hat', 'veg.hat', # 'ageM',                            # latent states
  'BetaA.S', 'BetaD.S', 'BetaV.S',                            # covariate effects
  'Mu.O', 'Epsilon.O', 'Sigma.O',                             # observation parameters
  'Gamma.S', 'Xi.S', 'Sigma.S',                               # random effects
  
  # Reproductive success model
  'Mu.B', 'Mu.Ri', 'Mu.Ra',                                   # mean reproductive success
  # 'BetaD.R', 'BetaV.R', 'BetaW.R',                          # covariate effects
  'EpsilonI.Ri', 'EpsilonT.Ri', 'EpsilonT.Ra', 'EpsilonT.B',  # random effects
  'SigmaI.Ri', 'SigmaT.Ri', 'SigmaT.Ra', 'SigmaT.B',          # random effects
  
  # Abundance model
  'ab', 'propF'
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

# start <- Sys.time()
# samples <- nimbleMCMC(code = myCode,
#                       data = myData,
#                       constants = myConst,
#                       inits = myInits,
#                       monitors = params,
#                       niter = niter,
#                       nburnin = nburnin,
#                       nchains = nchains,
#                       thin = nthin,
#                       samplesAsCodaMCMC = T,
#                       setSeed = seedMod)
# dur <- Sys.time() - start; dur
# beep(2)
# 
# # save MCMC output
# out.mcmc <- as.mcmc.list(samples)
# saveRDS(out.mcmc, 'results/IPM_CJSen_RSen_AB.rds', compress = 'xz')


## Run in parallel -------------------------------------------------------------

# function to run one chain inside cluster
runChain <- function(chainID, code, data, const, inits, params, niter, nburnin, nthin, seed){
  
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
                              "params", "nthin", "nburnin", "niter", "seedMod", "runChain"))

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

# combine & save
out.mcmc <- mcmc.list(samples)
saveRDS(out.mcmc, 'results/IPM_CJSen_RSen_AB.rds', compress = 'xz')


## Results ---------------------------------------------------------------------

library(coda)
library(MCMCvis)
library(corrplot)
library(ggplot2)
library(scales)

# # load results
# out.mcmc <- readRDS('results/IPM_CJS.rds')$out.mcmc
# summary(out.mcmc) # cannot handle NAs

# # find parameters generating NAs
# for(i in 1:ncol(out.mcmc[[1]])){
#   if(any(is.na(out.mcmc[[1]][,i]))){
#     message(paste0(colnames(out.mcmc[[1]])[i]))
#     print(out.mcmc[[1]][1:3,i])
#   }
# }

# summaries
MCMCsummary(out.mcmc, params = c('S'), n.eff = TRUE, round = 2)
MCMCsummary(out.mcmc, params = c('BetaA.S', 'BetaD.S', 'BetaV.S'), n.eff = TRUE, round = 2)
MCMCsummary(out.mcmc, params = c('Mu.O', 'Epsilon.O', 'Sigma.O'), n.eff = TRUE, round = 2)
MCMCsummary(out.mcmc, params = c('Sigma.S'), n.eff = TRUE, round = 2)

MCMCsummary(out.mcmc, params = c('Bt', 'Ra'), n.eff = TRUE, round = 2)
MCMCsummary(out.mcmc, params = c('BetaD.R', 'BetaV.R', 'BetaW.R'), n.eff = TRUE, round = 2)
MCMCsummary(out.mcmc, params = c('SigmaI.Ri', 'SigmaT.Ri', 'SigmaT.Ra', 'SigmaT.B'), n.eff = TRUE, round = 2)

MCMCsummary(out.mcmc, params = c('nYAF', 'nSA', 'nAD', 'nTOT'), n.eff = TRUE, round = 2)
MCMCsummary(out.mcmc, params = c('ab', 'propF'), n.eff = TRUE, round = 2)

# chainplots
MCMCtrace(out.mcmc, params = c('S'), pdf = FALSE)
MCMCtrace(out.mcmc, params = c('BetaA.S', 'BetaD.S', 'BetaV.S'), pdf = FALSE)
MCMCtrace(out.mcmc, params = c('Mu.O', 'Epsilon.O', 'Sigma.O'), pdf = FALSE)
MCMCtrace(out.mcmc, params = c('Sigma.S'), pdf = FALSE)

MCMCtrace(out.mcmc, params = c('Bt', 'Ra'), pdf = FALSE)
MCMCtrace(out.mcmc, params = c('BetaD.R', 'BetaV.R', 'BetaW.R'), pdf = FALSE)
MCMCtrace(out.mcmc, params = c('SigmaI.Ri', 'SigmaT.Ri', 'SigmaT.Ra', 'SigmaT.B'), pdf = FALSE)

MCMCtrace(out.mcmc, params = c('nYAF', 'nSA', 'nAD', 'nTOT'), pdf = FALSE)
MCMCtrace(out.mcmc, params = c('ab', 'propF'), pdf = FALSE)


## Compare model outputs -------------------------------------------------------

nYear <- myConst$nYear
nAge  <- myConst$nAge
nAgeC <- myConst$nAgeC

source('compareModels.R')
compareModels(nAge = nAge, nAgeC = nAgeC, nYear = nYear,
              postPaths = c("results/IPM_CJS_RS.rds", "results/IPM_CJS_RSen.rds"),
              modelNames = c("IPM/CJS/RS", "IPM/CJS/RSen"), plotFolder = c("figures/RSen"),
              returnSumData = TRUE)


## Extract parameter samples ---------------------------------------------------

source('extractParamSamples.R')
param.samples <- extractParamSamples(MCMCsamples = out.mcmc, saveList = TRUE)
# param.samples <- readRDS('results/paramSamples.rds')


## Plot population model -------------------------------------------------------

# posterior samples
out.mat <- as.matrix(out.mcmc)

# parameters to include
table.params <- c(
  paste0('nYAF[', 1:nYear, ']'),
  paste0('nSA[', 1:nYear, ']'),
  paste0('nAD[', rep(1:nAge, each = nYear), ', ', rep(1:nYear, times = nAge), ']'),
  paste0('nTOT[', 1:nYear, ']'))

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
# nYAF <- grep("^nYAF\\[", colnames(samples)); nYAF
# nSA <- grep("^nSA\\[", colnames(samples)); nSA
# nAD <- grep("^nAD\\[", colnames(samples)); nAD
# nTOT <- grep("^nTOT\\[", colnames(samples)); nTOT

nYAF <- grep("^nYAF\\[", colnames(out.mcmc[[1]])); nYAF
nSA <- grep("^nSA\\[", colnames(out.mcmc[[1]])); nSA
nAD <- grep("^nAD\\[", colnames(out.mcmc[[1]])); nAD
nTOT <- grep("^nTOT\\[", colnames(out.mcmc[[1]])); nTOT

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
df <- expand.grid(age = 1:22, year = 1:16)
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
  geom_ribbon(aes(ymin = sLCI, ymax = sUCI, fill = ageC), alpha = 0.2) +
  geom_line(aes(colour = ageC), linewidth = 1, show.legend = F) +
  labs(x = "Year", y = "Survival", fill = "Age class") +
  # scale_x_continuous(breaks = pretty_breaks()) +
  scale_y_continuous(breaks = pretty_breaks()) +
  theme_bw()

# ggsave("figures/IPM_CJS.jpeg", scale = 1, width = 18.0, height = 9.0, units = c("cm"), dpi = 600)

