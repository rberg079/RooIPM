# 8 April 2025
# Process model for roo IPM
# 1 census per year, on Sept 1 ish

## Set up ----------------------------------------------------------------------

# set toggles
testRun <- TRUE
use_dCJS <- TRUE
parallelRun <- TRUE
envEffectsS <- TRUE
envEffectsR <- TRUE
ageClasses <- 12
splitCovs <- 2
splitREs <- 1

# load packages
library(tidyverse)
library(lubridate)
library(beepr)
library(boot)
library(coda)
library(here)
library(nimble)
library(nimbleEcology)
library(parallel)
library(corrplot)
library(MCMCvis)

# load data
source('R/wrangleData_en.R')
enData <- wrangleData_en(
  dens.data = "data/WPNP_Methods_Results_January2026.xlsx",
  veg.data  = "data/biomass data April 2009 - July 2025_updated Feb2026.xlsx",
  wea.data  = "data/Prom_Weather_2008-2023_updated Jan2026 RB.xlsx",
  wind.data = "data/POWER_Point_Daily_20080101_20260331_10M.csv",
  obs.data  = "data/PromObs_2008-2024.xlsx",
  list.data = "data/PromlistAllNov25.xlsx")

source('R/wrangleData_sv.R')
svData <- wrangleData_sv(
  surv.data = "data/PromSurvivalNov25_RB.xlsx",
  yafs.data = "data/RSmainRB_May26.xlsx",
  ageClasses = ageClasses, known.age = TRUE, from2012 = FALSE, splitCovs = splitCovs)

source('R/wrangleData_rs.R')
rsData <- wrangleData_rs(
  rs.data = "data/RSmainRB_May26.xlsx",
  ageClasses = ageClasses, known.age = TRUE, cum.surv = FALSE)

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
                vegE = enData$vegE)

myConst <- list(nYear = svData$nYear,
                nAge = rsData$nAge+1,
                
                nID.S = svData$nID,
                nAgeC.S = svData$nAgeC.S,
                age.S = svData$age.S,
                ageC.S = svData$ageC.S,
                idx.sv = svData$idx.sv,
                
                nB = rsData$nB,
                nR = rsData$nR,
                nID.R = rsData$nID,
                id.R = rsData$id.R,
                year.B = rsData$year.B,
                year.R = rsData$year.R,
                age.B = rsData$age.B,
                age.R = rsData$age.R,
                ageC.R = rsData$ageC.R,
                nAgeC.R = rsData$nAgeC.R,
                
                first = svData$first,
                last = svData$last,
                
                densM = enData$densM,
                noVeg = enData$noVeg,
                noProp = enData$noProp,
                nNoVeg = enData$nNoVeg,
                nNoProp = enData$nNoProp,
                
                use_dCJS = use_dCJS,
                envEffectsS = envEffectsS,
                envEffectsR = envEffectsR,
                ageClasses = ageClasses,
                splitCovs = splitCovs,
                splitREs = splitREs)

# conditionally add dummy variables
if(splitCovs == 3){
  myConst <- c(myConst, list(dummyY = svData$dummyY, dummyP = svData$dummyP, dummyO = svData$dummyO))
}else if(splitCovs == 2){
  myConst <- c(myConst, list(dummyY = svData$dummyY, dummyO = svData$dummyO))
}else if(splitCovs == 1){
  myConst <- c(myConst, list(dummy = svData$dummy))
}

## Assemble --------------------------------------------------------------------

source('R/writeCode.R')
myCode <- writeCode()

nchains   <- 8
seedMod   <- c(30, 31, 32, 33, 34, 35, 36, 37)
seedInits <- 38

# assign initial values
source('R/simulateInits.R')
set.seed(seedInits)
myInits <- list()
for(c in 1:nchains){
  myInits[[c]] <- simulateInits(
    dens = myData$dens,
    veg = myData$veg,
    propF = myData$propF,
    knownStates = svData$state,
    nYear = myConst$nYear,
    nAge = myConst$nAge,
    nB = myConst$nB,
    nR = myConst$nR,
    nID.R = myConst$nID.R,
    ageClasses = ageClasses,
    year.B = myConst$year.B,
    year.R = myConst$year.R,
    id.R = myConst$id.R,
    age.B = myConst$age.B,
    age.R = myConst$age.R,
    ageC.R = myConst$ageC.R,
    ageC.S = myConst$ageC.S,
    envEffectsS = envEffectsS,
    envEffectsR = envEffectsR,
    splitCovs = splitCovs,
    splitREs = splitREs
    )
}

# select parameters to monitors
params <- c(
  # Population model
  'S', 'BR', 'sPY', 'sYF', 'sSA', 'sAD',
  'nYF', 'nSA', 'nAD', 'nTOT',
  
  # Survival model
  'Mu.S',
  'Mu.O', 'EpsilonT.O', 'SigmaT.O',
  
  # Reproductive success model
  'Mu.B', 'Mu.R', 
  'EpsilonT.B', 'EpsilonI.R', 'EpsilonT.R', 
  'SigmaT.B', 'SigmaI.R', 'SigmaT.R', 
  
  # Density model
  'propF'
)

# conditionally add covariate effects
if(envEffectsS){
  if(splitCovs == 3){
    params <- c(params, 'BetaD.Sy', 'BetaD.Sp', 'BetaD.So', 'BetaV.Sy', 'BetaV.Sp', 'BetaV.So')
  }else if(splitCovs == 2){
    params <- c(params, 'BetaD.Sy', 'BetaD.So', 'BetaV.Sy', 'BetaV.So')
  }else if(splitCovs == 1){
    params <- c(params, 'BetaD.S', 'BetaV.S')
  }
} 

if(envEffectsR){params <- c(params, 'BetaD.R')}
if(envEffectsS || envEffectsR){params <- c(params, 'dens.true', 'veg.true')}

# conditionally add random effects
if(splitREs == 3){
    params <- c(params, 'EpsilonT.Sy', 'EpsilonT.Sp', 'EpsilonT.So', 'SigmaT.Sy', 'SigmaT.Sp', 'SigmaT.So')
  }else if(splitREs == 2){
    params <- c(params, 'EpsilonT.Sy', 'EpsilonT.So', 'SigmaT.Sy', 'SigmaT.So')
  }else if(splitREs == 1){
    params <- c(params, 'EpsilonT.S', 'SigmaT.S')
  }

# select MCMC settings
if(testRun){
  nthin   <- 1
  nburnin <- 0
  niter   <- 10
}else{
  nthin   <- 20
  nburnin <- 60000
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
saveRDS(out.mcmc, 'results/IPM_CJSen_RSen_AB_DynDens_dCJS_12_noW_25BR_Dave2Covs_stochV_8chains.rds', compress = 'xz')


## Results ---------------------------------------------------------------------

# # load results
# out.mcmc <- readRDS('results/IPM_CJSen_RSen_AB_DynDens_dCJS_12_noW_stochV_25&BR_Dave3CovsII.rds')
# summary(out.mcmc) # cannot handle NAs

# # find parameters generating NAs
# for(i in 1:ncol(out.mcmc[[1]])){
#   if(any(is.na(out.mcmc[[1]][,i]))){
#     message(paste0(colnames(out.mcmc[[1]])[i]))
#     print(out.mcmc[[1]][1:3,i])
#   }
# }

# ## summaries
# # population model
# MCMCsummary(out.mcmc, params = c('S'), n.eff = TRUE, round = 2)
# MCMCsummary(out.mcmc, params = c('BR', 'sPY', 'sYF', 'sSA', 'sAD'), n.eff = TRUE, round = 2)
# MCMCsummary(out.mcmc, params = c('nYF', 'nSA', 'nAD', 'nTOT', 'propF'), n.eff = TRUE, round = 2)
# 
# # survival model
# MCMCsummary(out.mcmc, params = c('Mu.S', 'Mu.O', 'EpsilonT.O', 'SigmaT.O'), n.eff = TRUE, round = 2)
# 
# if(envEffectsS){
#   if(splitCovs == 3){
#     MCMCsummary(out.mcmc, params = c('BetaD.Sy', 'BetaD.So', 'BetaD.Sp', 'BetaV.Sy', 'BetaV.So', 'BetaV.Sp'), n.eff = TRUE, round = 2, pg0 = TRUE)
#   } else if(splitCovs == 2){
#     MCMCsummary(out.mcmc, params = c('BetaD.Sy', 'BetaD.So', 'BetaV.Sy', 'BetaV.So'), n.eff = TRUE, round = 2, pg0 = TRUE)
#   } else if(splitCovs == 1){
#     MCMCsummary(out.mcmc, params = c('BetaD.S', 'BetaV.S'), n.eff = TRUE, round = 2, pg0 = TRUE)
#   }
# }
# 
# if(splitREs == 3){
#   MCMCsummary(out.mcmc, params = c('EpsilonT.Sy', 'EpsilonT.Sp', 'EpsilonT.So', 'SigmaT.Sy', 'SigmaT.Sp', 'SigmaT.So'), n.eff = TRUE, round = 2, pg0 = TRUE)
# }else if(splitREs == 2){
#   MCMCsummary(out.mcmc, params = c('EpsilonT.Sy', 'EpsilonT.So', 'SigmaT.Sy', 'SigmaT.So'), n.eff = TRUE, round = 2, pg0 = TRUE)
# }else if(splitREs == 1){
#   MCMCsummary(out.mcmc, params = c('EpsilonT.S', 'SigmaT.S'), n.eff = TRUE, round = 2, pg0 = TRUE)
# }
# 
# # reproductive success model
# MCMCsummary(out.mcmc, params = c('Mu.B', 'Mu.R'), n.eff = TRUE, round = 2)
# MCMCsummary(out.mcmc, params = c('EpsilonT.B', 'EpsilonI.R', 'EpsilonT.R'), n.eff = TRUE, round = 2)
# MCMCsummary(out.mcmc, params = c('SigmaT.B', 'SigmaI.R', 'SigmaT.R'), n.eff = TRUE, round = 2)
# 
# if(envEffectsR){
#   MCMCsummary(out.mcmc, params = c('BetaD.R'), n.eff = TRUE, round = 2, pg0 = TRUE)
# }
# 
# # latent true environment
# if(envEffectsS || envEffectsR){
#   MCMCsummary(out.mcmc, params = c('dens.true', 'veg.true'), n.eff = TRUE, round = 2)
# }
# 
# ## chainplots
# # population model
# MCMCtrace(out.mcmc, params = c('S'), pdf = FALSE)
# MCMCtrace(out.mcmc, params = c('BR', 'sPY', 'sYF', 'sSA', 'sAD'), pdf = FALSE)
# MCMCtrace(out.mcmc, params = c('nYF', 'nSA', 'nAD', 'nTOT', 'propF'), pdf = FALSE)
# 
# # survival model
# MCMCtrace(out.mcmc, params = c('Mu.S', 'EpsilonT.S', 'SigmaT.S'), pdf = FALSE)
# MCMCtrace(out.mcmc, params = c('Mu.O', 'EpsilonT.O', 'SigmaT.O'), pdf = FALSE)
# 
# if(envEffectsS){
#   if(splitCovs == 3){
#     MCMCtrace(out.mcmc, params = c('BetaD.Sy', 'BetaD.So', 'BetaD.Sp', 'BetaV.Sy', 'BetaV.So', 'BetaV.Sp'), pdf = FALSE)
#   }else if(splitCovs == 2){
#     MCMCtrace(out.mcmc, params = c('BetaD.Sy', 'BetaD.So', 'BetaV.Sy', 'BetaV.So'), pdf = FALSE)
#   }else if(splitCovs == 1){
#     MCMCtrace(out.mcmc, params = c('BetaD.S', 'BetaV.S'), pdf = FALSE)
#   }
# }
# 
# if(splitREs == 3){
#   MCMCtrace(out.mcmc, params = c('EpsilonT.Sy', 'EpsilonT.Sp', 'EpsilonT.So', 'SigmaT.Sy', 'SigmaT.Sp', 'SigmaT.So'), pdf = FALSE)
# }else if(splitREs == 2){
#   MCMCtrace(out.mcmc, params = c('EpsilonT.Sy', 'EpsilonT.So', 'SigmaT.Sy', 'SigmaT.So'), pdf = FALSE)
# }else if(splitREs == 1){
#   MCMCtrace(out.mcmc, params = c('EpsilonT.S', 'SigmaT.S'), pdf = FALSE)
# }
# 
# # reproductive success model
# MCMCtrace(out.mcmc, params = c('Mu.B', 'Mu.R'), pdf = FALSE)
# MCMCtrace(out.mcmc, params = c('EpsilonT.B', 'EpsilonI.R', 'EpsilonT.R'), pdf = FALSE)
# MCMCtrace(out.mcmc, params = c('SigmaT.B', 'SigmaI.R', 'SigmaT.R'), pdf = FALSE)
# 
# if(envEffectsR){
#   MCMCtrace(out.mcmc, params = c('BetaD.R'), pdf = FALSE)
# }
# 
# # latent true environment
# if(envEffectsS || envEffectsR){
#   MCMCtrace(out.mcmc, params = c('dens.true', 'veg.true'), pdf = FALSE)
# }
# 
# ## export
# MCMCtrace(out.mcmc, Rhat = TRUE, pdf = TRUE, filename = 'figures/MCMCtrace.pdf')


## Compare model outputs -------------------------------------------------------

# nYear   <- myConst$nYear
# nAgeC.S <- myConst$nAgeC.S
# 
# source('R/compareModels.R')
# compareModels(nYear = nYear,
#               nAgeC.S = nAgeC.S,
#               postPaths = c(
#                 "results/IPM_CJSen_RSen_AB_DynDens_dCJS_12_noW_25BR_Dave2Covs_stochV_8chains.rds",
#                 "results/IPM_CJSen_RSen_AB_DynDens_dCJS_12_noW_25BR_Dave3Covs3REs_stochV_8chains.rds"
#               ),
#               modelNames = c(
#                 "IPM_Dave2Covs",
#                 "IPM_Dave3Covs3REs"
#               ),
#               plotFolder = c("figures/densityChecks/final2"),
#               returnSumData = TRUE)


## Extract parameter samples ---------------------------------------------------

# source('R/extractParamSamples.R')
# out.mcmc <- readRDS('results/IPM_CJSen_RSen_AB_DynDens_dCJS_12_noW_25BR_Dave2Covs_stochV_8chains.rds')
# paramSamples <- extractParamSamples(MCMCsamples = out.mcmc, saveList = TRUE)


## Calculate sensitivities & elasticities --------------------------------------

# source('R/calculateSensitivities.R')
# sensitivities <- calculateSensitivities(paramSamples = paramSamples)


## Run random design transient LTRE --------------------------------------------

# source('R/runLTRE_random.R')
# LTREresults <- runLTRE_random(paramSamples = paramSamples, sensitivities = sensitivities)

