# 8 April 2025
# Process model for roo IPM
# 1 census per year, on Sept 1 ish

## Set up ----------------------------------------------------------------------

# set toggles
testRun <- FALSE
parallelRun <- TRUE
envEffectsS <- TRUE
envEffectsR <- TRUE
ageClasses <- 12
use_dCJS <- TRUE

# load packages
library(tidyverse)
library(lubridate)
library(beepr)
library(here)
library(boot)
library(coda)
library(nimble)
library(nimbleEcology)
library(parallel)

# load data
source('wrangleData_en.R')
enData <- wrangleData_en(dens.data = "data/abundanceData_Proteus.csv",
                         veg.data  = "data/biomass data April 2009 - Jan 2025_updated Feb2025.xlsx",
                         wea.data  = "data/Prom_Weather_2008-2023_updated Jan2025 RB.xlsx",
                         wind.data = "data/POWER_Point_Daily_20080101_20241231_10M.csv",
                         obs.data  = "data/PromObs_2008-2024.xlsx",
                         list      = "data/PromlistAllOct24.xlsx")

source('wrangleData_sv.R')
svData <- wrangleData_sv(surv.data = "data/PromSurvivalOct24.xlsx",
                         yafs.data = "data/RSmainRB_Mar25.xlsx",
                         ageClasses = ageClasses, known.age = TRUE)

source('wrangleData_rs.R')
rsData <- wrangleData_rs(rs.data = "data/RSmainRB_Mar25.xlsx",
                         obs.data = "data/PromObs_2008-2024.xlsx",
                         ageClasses = ageClasses, known.age = TRUE, cum.surv = FALSE)

# NAs in age.S before first capture were throwing an error at model defining step!
# replacing NAs with dummy integer 1 seems to have solved it (to move to wrangling)
svData$age.S[is.na(svData$age.S)] <- 1

# create Nimble lists
myData  <- list(obs = svData$obs,
                state = svData$state,
                
                B = rsData$B,
                R = rsData$R,
                
                area = enData$area,
                propF = enData$propF,
                dens = enData$dens,
                densE = enData$densE,
                veg = enData$veg,
                vegE = enData$vegE,
                win = enData$win)

myConst <- list(nYear = svData$nYear,
                nAge = rsData$nAge+1,
                
                nID.S = svData$nID,
                nID.S.switch = min(which(svData$first == svData$nYear - 1)),
                nAgeC.S = svData$nAgeC.S,
                age.S = svData$age.S,
                ageC.S = svData$ageC.S,
                
                nB = rsData$nB,
                nR = rsData$nR,
                nID.R = rsData$nID,
                id.R = rsData$id.R,
                year.B = rsData$year.B,
                year.R = rsData$year.R,
                age.R = rsData$age.R,
                ageC.R = rsData$ageC.R,
                nAgeC.R = rsData$nAgeC.R,
                
                dummy = svData$dummy,
                first = svData$first,
                last = svData$last,
                
                densM = enData$densM,
                noVeg = enData$noVeg,
                noWin = enData$noWin,
                noProp = enData$noProp,
                nNoVeg = enData$nNoVeg,
                nNoWin = enData$nNoWin,
                nNoProp = enData$nNoProp,
                
                envEffectsS = envEffectsS,
                envEffectsR = envEffectsR,
                ageClasses = ageClasses,
                use_dCJS = use_dCJS)


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
    dens = myData$dens,
    veg = myData$veg,
    win = myData$win,
    propF = myData$propF,
    knownStates = svData$state,
    nYear = myConst$nYear,
    nAge = myConst$nAge,
    nR = myConst$nR,
    nID.R = myConst$nID.R,
    ageClasses = ageClasses,
    year.R = myConst$year.R,
    id.R = myConst$id.R,
    age.R = myConst$age.R,
    ageC.R = myConst$ageC.R,
    ageC.S = myConst$ageC.S,
    envEffectsR = TRUE,
    envEffectsS = TRUE
    )
}

# select parameters to monitors
params <- c(
  # Population model
  'S', 'Bt', 'sPY', 'sYF', 'sSA', 'sAD',
  'nYF', 'nSA', 'nAD', 'nTOT',
  
  # Survival model
  'Mu.S', 'EpsilonT.S', 'SigmaT.S',
  'Mu.O', 'EpsilonT.O', 'SigmaT.O',
  
  # Reproductive success model
  'Mu.B', 'Mu.R', 
  'EpsilonT.B', 'EpsilonI.R', 'EpsilonT.R', 
  'SigmaT.B', 'SigmaI.R', 'SigmaT.R', 
  
  # Abundance model
  'propF'
)

# conditionally add covariate effects
if(envEffectsS){params <- c(params, 'BetaD.S', 'BetaV.S', 'BetaW.S')}
if(envEffectsR){params <- c(params, 'BetaD.R')}
if(envEffectsS || envEffectsR){params <- c(params, 'dens.true', 'veg.true', 'win.true')}

# select MCMC settings
if(testRun){
  nthin   <- 1
  nburnin <- 0
  niter   <- 10
}else{
  nthin   <- 20
  nburnin <- 20000
  niter   <- nburnin + 1000*nthin
}


## Run model -------------------------------------------------------------------

if(parallelRun){
  # function to run one chain inside cluster
  runChain <- function(chainID, code, data, const, inits, params,
                       niter, nburnin, nthin, seed){
    
    library(nimble)
    library(nimbleEcology)
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
saveRDS(out.mcmc, 'results/IPM_CJSen_RSen_AB_DynDens_dCJS_12.rds', compress = 'xz')


## Results ---------------------------------------------------------------------

library(coda)
library(MCMCvis)
library(corrplot)
library(ggplot2)
library(scales)

# # load results
# out.mcmc <- readRDS('results/IPM_CJSen_RSen_AB_DynDens_dCJS_12.rds')
# summary(out.mcmc) # cannot handle NAs

# # find parameters generating NAs
# for(i in 1:ncol(out.mcmc[[1]])){
#   if(any(is.na(out.mcmc[[1]][,i]))){
#     message(paste0(colnames(out.mcmc[[1]])[i]))
#     print(out.mcmc[[1]][1:3,i])
#   }
# }

# # summaries
# MCMCsummary(out.mcmc, params = c('S'), n.eff = TRUE, round = 2)
# MCMCsummary(out.mcmc, params = c('Mu.S', 'EpsilonT.S', 'SigmaT.S'), n.eff = TRUE, round = 2)
# if(envEffectsS){MCMCsummary(out.mcmc, params = c('BetaD.S', 'BetaV.S', 'BetaW.S'), n.eff = TRUE, round = 2, pg0 = T)}
# MCMCsummary(out.mcmc, params = c('Mu.O', 'EpsilonT.O', 'SigmaT.O'), n.eff = TRUE, round = 2)
# 
# MCMCsummary(out.mcmc, params = c('Mu.R', 'Bt', 'sPY'), n.eff = TRUE, round = 2)
# if(envEffectsR){MCMCsummary(out.mcmc, params = c('BetaD.R'), n.eff = TRUE, round = 2, pg0 = T)}
# MCMCsummary(out.mcmc, params = c('SigmaI.R', 'SigmaT.R', 'SigmaT.B'), n.eff = TRUE, round = 2)
# 
# MCMCsummary(out.mcmc, params = c('nYF', 'nSA', 'nAD', 'nTOT', 'propF'), n.eff = TRUE, round = 2)
# 
# # chainplots
# MCMCtrace(out.mcmc, params = c('S'), pdf = FALSE)
# MCMCsummary(out.mcmc, params = c('Mu.S', 'EpsilonT.S', 'SigmaT.S'), n.eff = TRUE, round = 2)
# if(envEffectsS){MCMCtrace(out.mcmc, params = c('BetaD.S', 'BetaV.S', 'BetaW.S'), pdf = FALSE)}
# MCMCtrace(out.mcmc, params = c('Mu.O', 'EpsilonT.O', 'SigmaT.O'), pdf = FALSE)
# 
# MCMCtrace(out.mcmc, params = c('Bt', 'sPY'), pdf = FALSE)
# if(envEffectsR){MCMCtrace(out.mcmc, params = c('BetaD.R'), pdf = FALSE)}
# MCMCtrace(out.mcmc, params = c('SigmaI.R', 'SigmaT.R', 'SigmaT.B'), pdf = FALSE)
# 
# MCMCtrace(out.mcmc, params = c('nYF', 'nSA', 'nAD', 'nTOT', 'propF'), pdf = FALSE)
# 
# MCMCtrace(out.mcmc, Rhat = T, pdf = T, filename = 'results/MCMCtrace.pdf')


## Compare model outputs -------------------------------------------------------

# nYear   <- myConst$nYear
# nAgeC.S <- myConst$nAgeC.S
# 
# source('compareModels.R')
# compareModels(nYear = nYear,
#               nAgeC.S = nAgeC.S,
#               postPaths = c(
#                 "results/IPM_CJSen_RSen_AB_DynDens_base.rds",
#                 "results/IPM_CJSen_RSen_AB_DynDens_dCJS.rds",
#                 "results/IPM_CJSen_RSen_AB_DynDens_dCJS_12.rds",
#                 "results/IPM_CJSen_RSen_AB_DynDens_dCJS_20.rds"
#               ),
#               modelNames = c(
#                 "base",
#                 "dCJS_6",
#                 "dCJS_12",
#                 "dCJS_20"
#               ),
#               plotFolder = c("figures/ageClasses"),
#               returnSumData = TRUE)


## Extract parameter samples ---------------------------------------------------

# source('extractParamSamples.R')
# paramSamples <- extractParamSamples(MCMCsamples = out.mcmc, saveList = TRUE)
# # paramSamples <- readRDS('results/paramSamples.rds')

