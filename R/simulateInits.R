#' Simulate initial values for IPM
#'
#' @param dens vector of length nYear of annual population density.
#' @param veg vector of length nYear of annual available vegetation.
#' @param propF vector of length nYear of proportion of kangaroos observations that are of females.
#' @param knownStates matrix of same dimensions as capture history of information on known states.
#' @param nYear integer. Number of time steps, or years, in the model. nYear = 18 by default.
#' @param nAge integer. Number of ages, or maximum age, in the model. nAge = 19 by default.
#' @param nB integer. Number of reproductive events in the birth rate model.
#' @param nR integer. Number of reproductive events in the reproductive success model.
#' @param nID.R integer. Number of unique female kangaroos in the reproductive success model.
#' @param ageClasses integer. Number of age classes to be considered. ageClasses = 20 by default.
#' @param year.B integer vector. Year of each reproductive event in the birth rate analysis.
#' @param year.R integer vector. Year of each reproductive event in the reproductive success analysis.
#' @param id.R integer vector. Maternal ID of each reproductive event in the reproductive success analysis.
#' @param age.B integer vector. Maternal age of each reproductive event in the birth rate analysis.
#' @param age.R integer vector. Maternal age of each reproductive event in the reproductive success analysis.
#' @param ageC.R integer vector. Vector mapping age classes to actual ages in the reproductive success model.
#' @param envEffectsS logical. If TRUE, environmental covariates are included in CJS model.
#' @param envEffectsR logical. If TRUE, environmental covariates are included in RS model.
#' @param splitCovs integer. Number of age-group-specific covariate effects to fit in the survival model. splitCovs = 1 by default.
#' @param splitREs integer. Number of age-group-specific random effects of year to fit in the survival model. splitREs = 1 by default.
#'
#' @returns a list containing all initial values needed for the IPM.
#' @export
#'
#' @examples

simulateInits <- function(dens, veg, propF, knownStates, 
                          nYear = 18, nAge = 19, nB, nR, nID.R, ageClasses = 20,
                          year.B, year.R, id.R, age.B, age.R, ageC.R, ageC.S,
                          envEffectsS = TRUE, envEffectsR = TRUE, splitCovs = 1, splitREs = 1){
  
  # # for testing purposes
  # library(readxl)
  # library(tidyverse)
  # 
  # ageClasses <- 12
  # source('R/wrangleData_en.R')
  # enData <- wrangleData_en(
  #   dens.data = "data/WPNP_Methods_Results_January2026.xlsx",
  #   veg.data  = "data/biomass data April 2009 - July 2025_updated Feb2026.xlsx",
  #   wea.data  = "data/Prom_Weather_2008-2023_updated Jan2026 RB.xlsx",
  #   wind.data = "data/POWER_Point_Daily_20080101_20260331_10M.csv",
  #   obs.data  = "data/PromObs_2008-2024.xlsx",
  #   list.data = "data/PromlistAllNov25.xlsx")
  # 
  # source('R/wrangleData_sv.R')
  # svData <- wrangleData_sv(
  #   surv.data = "data/PromSurvivalNov25_RB.xlsx",
  #   yafs.data = "data/RSmainRB_May26.xlsx",
  #   ageClasses = ageClasses, known.age = TRUE, from2012 = FALSE, splitCovs = 2)
  # 
  # source('R/wrangleData_rs.R')
  # rsData <- wrangleData_rs(
  #   rs.data = "data/RSmainRB_May26.xlsx",
  #   ageClasses = ageClasses, known.age = TRUE, cum.surv = FALSE)
  # 
  # dens <- enData$dens
  # veg <- enData$veg
  # propF <- enData$propF
  # knownStates <- svData$state
  # 
  # nYear <- 18
  # nAge <- 19
  # 
  # nB <- rsData$nB
  # nR <- rsData$nR
  # nID.R <- rsData$nID.R
  # year.B <- rsData$year.B
  # year.R <- rsData$year.R
  # id.R <- rsData$id.R
  # age.B <- rsData$age.B
  # age.R <- rsData$age.R
  # ageC.R <- rsData$ageC.R
  # ageC.S <- svData$ageC.S
  # 
  # envEffectsR <- TRUE
  # envEffectsS <- TRUE
  # splitCovs <- 2
  # splitREs <- 1
  
  
  ## Simulate latent states for input data -------------------------------------
  
  ## Survival model
  # missing values
  noDens <- is.na(dens)
  noVeg  <- is.na(veg)
  noProp <- is.na(propF)
  
  # nNoDens <- length(which(noDens))
  nNoVeg  <- length(which(noVeg))
  nNoProp <- length(which(noProp))
  
  dens  <- round(ifelse(noDens, rnorm(1, 3.0, 1), dens), 2)
  veg   <- round(ifelse(noVeg, rnorm(1, 0, .1), veg), 4)
  propF <- round(ifelse(noProp, pmax(pmin(rnorm(1, .7, .1), 0.99), 0.4), propF), 4)
  
  # true environment
  dens.true <- dens
  dens.cov  <- dens - mean(dens)
  veg.true  <- veg
  
  # latent states
  state <- knownStates
  
  for(i in 1:nrow(state)){
    first <- min(which(!is.na(state[i, ]))) # first capture
    state[i, is.na(state[i, ])] <- 1        # 1 for any unknown states
    state[i, 1:(first-1)] <- 0              # 0 for any NA prior to first
  }
  
  ## Simulate vital rate covariate effects -------------------------------------
  
  nAgeC.S <- max(ageC.S)
  nAgeC.R <- max(ageC.R)
  
  ## Survival model
  if(envEffectsS){
    if(splitCovs == 3){
      BetaD.Sy <- runif(1, -1, 1)
      BetaD.Sp <- runif(1, -1, 1)
      BetaD.So <- runif(1, -1, 1)
      BetaV.Sy <- runif(1, -1, 1)
      BetaV.Sp <- runif(1, -1, 1)
      BetaV.So <- runif(1, -1, 1)
    }else if(splitCovs == 2){
      BetaD.Sy <- runif(1, -1, 1)
      BetaD.So <- runif(1, -1, 1)
      BetaV.Sy <- runif(1, -1, 1)
      BetaV.So <- runif(1, -1, 1)
    }else if(splitCovs == 1){
      BetaD.S <- runif(1, -1, 1)
      BetaV.S <- runif(1, -1, 1)
    }
  }
  
  ## Reproductive success model
  if(envEffectsR){
    BetaD.R <- runif(1, -1, 1)
  }
  
  # dummy variable
  # to target covariate effects
  if(splitCovs > 1){
    if(ageClasses == 6){
      dummyY = c(1, rep(0, 5))
      dummyP = c(0, rep(1, 4), 0)
      dummyO = c(rep(0, 5), 1)
    }else if(ageClasses == 12){
      dummyY = c(1, rep(0, 12))
      dummyP = c(0, rep(1, 9), rep(0, 3))
      dummyO = c(rep(0, 10), rep(1, 3))
    }else if(ageClasses == 20){
      dummyY = c(1, rep(0, 19))
      dummyP = c(0, rep(1, 9), rep(0, 10))
      dummyO = c(rep(0, 10), rep(1, 10))
    }
    
  }else if(splitCovs == 1){
    if(ageClasses == 6){
      dummy = c(1, rep(0,4), 1)
    }else if(ageClasses == 12){
      dummy = c(1, rep(0,8), rep(1,4))
    }else if(ageClasses == 20){
      dummy = c(1, rep(0,8), rep(1,11))
    }
  }
  
  
  ## Simulate vital rate random effects ----------------------------------------
  
  ## Survival model
  if(splitREs == 3){
    XiT.Sy <- rnorm(nYear-1, 0, 1)
    SigmaT.Sy <- runif(1, .5, 2)
    EpsilonT.Sy <- XiT.Sy * SigmaT.Sy
    
    XiT.Sp <- rnorm(nYear-1, 0, 1)
    SigmaT.Sp <- runif(1, .5, 2)
    EpsilonT.Sp <- XiT.Sp * SigmaT.Sp
    
    XiT.So <- rnorm(nYear-1, 0, 1)
    SigmaT.So <- runif(1, .5, 2)
    EpsilonT.So <- XiT.So * SigmaT.So
    
  }else if(splitREs == 2){
    XiT.Sy <- rnorm(nYear-1, 0, 1)
    SigmaT.Sy <- runif(1, .5, 2)
    EpsilonT.Sy <- XiT.Sy * SigmaT.Sy
    
    XiT.So <- rnorm(nYear-1, 0, 1)
    SigmaT.So <- runif(1, .5, 2)
    EpsilonT.So <- XiT.So * SigmaT.So
    
  }else if(splitREs == 1){
    XiT.S <- rnorm(nYear-1, 0, 1)
    SigmaT.S <- runif(1, .5, 2)
    EpsilonT.S <- XiT.S * SigmaT.S
  }
  
  ## Reproductive success model
  XiI.R <- rnorm(nID.R, 0, 1)
  XiT.R <- rnorm(nYear-1, 0, 1)
  XiT.B <- rnorm(nYear-1, 0, 1)

  SigmaI.R <- runif(1, .5, 2)
  SigmaT.R <- runif(1, .5, 2)
  SigmaT.B <- runif(1, .5, 2)
  
  EpsilonI.R <- XiI.R * SigmaI.R
  EpsilonT.R <- XiT.R * SigmaT.R
  EpsilonT.B <- XiT.B * SigmaT.B

  
  ## Simulate yearly vital rates -----------------------------------------------
  
  ## Survival model
  Mu.S <- c(rep(runif(nAgeC.S, 0.1, 0.9)))
  S <- matrix(NA, nrow = nAgeC.S, ncol = nYear-1)
  
  for(a in 1:nAgeC.S){
    for(t in 1:(nYear - 1)){
      
      # intercepts
      logit.S <- qlogis(Mu.S[a])
      
      # covariate effects
      if(envEffectsS){
        if(splitCovs == 3){
          logit.S <- logit.S + 
            BetaD.Sy * dens.cov[t] * dummyY[a] + 
            BetaD.Sp * dens.cov[t] * dummyP[a] + 
            BetaD.So * dens.cov[t] * dummyO[a] +
            BetaV.Sy * veg.true[t] * dummyY[a] + 
            BetaV.Sp * veg.true[t] * dummyP[a] + 
            BetaV.So * veg.true[t] * dummyO[a]
        }else if(splitCovs == 2){
          logit.S <- logit.S + 
            BetaD.Sy * dens.cov[t] * dummyY[a] + 
            BetaD.So * dens.cov[t] * dummyO[a] +
            BetaV.Sy * veg.true[t] * dummyY[a] + 
            BetaV.So * veg.true[t] * dummyO[a]
        }else if(splitCovs == 1){
          logit.S <- logit.S + 
            BetaD.S * dens.cov[t] * dummy[a] + 
            BetaV.S * veg.true[t] * dummy[a]
        }
      }
      
      # random effects
      if(splitREs == 3){
        logit.S <- logit.S + 
          EpsilonT.Sy[t] * dummyY[a] + 
          EpsilonT.Sp[t] * dummyP[a] + 
          EpsilonT.So[t] * dummyO[a]
      }else if(splitREs == 2){
        logit.S <- logit.S + 
          EpsilonT.Sy[t] * dummyY[a] + 
          EpsilonT.So[t] * dummyO[a]
      }else if(splitREs == 1){
        logit.S <- logit.S +
          EpsilonT.S[t]
      }
      
      S[a, t] <- plogis(logit.S)
    }
  }
  
  ## Reproductive success model
  Mu.B <- c(rep(runif(nAgeC.R, 0.1, 0.9)))

  Bi <- numeric(nB)
  for(x in 1:nB){
    Bi[x] <- plogis(
      qlogis(Mu.B[ageC.R[age.B[x]]]) + 
        EpsilonT.B[year.B[x]])
  }
  
  Ba <- matrix(0, nrow = nAgeC.R, ncol = nYear-1)
  for(a in 1:nAgeC.R){
    for(t in 1:(nYear-1)){
      Ba[a, t] <- plogis(
        qlogis(Mu.B[a]) + 
          EpsilonT.B[t])
    }
  }
  
  # age-specific baseline reproductive success
  Mu.R <- c(rep(runif(nAgeC.R, 0.1, 0.9)))
  
  Ri <- numeric(nR)
  for(x in 1:nR){
    if(envEffectsR){
      Ri[x] <- plogis(
        qlogis(Mu.R[ageC.R[age.R[x]]]) +
          BetaD.R * dens.cov[year.R[x]] +
          EpsilonI.R[id.R[x]] +
          EpsilonT.R[year.R[x]])
    }else{
      Ri[x] <- plogis(
        qlogis(Mu.R[ageC.R[age.R[x]]]) +
          EpsilonI.R[id.R[x]] +
          EpsilonT.R[year.R[x]])
    }
  }
  
  Ra <- matrix(0, nrow = nAgeC.R, ncol = nYear-1)
  for(a in 1:nAgeC.R) {
    for(t in 1:(nYear-1)) {
      if(envEffectsR){
        Ra[a, t] <- plogis(
          qlogis(Mu.R[a]) +
            BetaD.R * dens.cov[t] +
            EpsilonT.R[t])
      }else{
        Ra[a, t] <- plogis(
          qlogis(Mu.R[a]) +
            EpsilonT.R[t])
      }
    }
  }
  
  ## Population model
  # priors for survival
  sYF <- S[1, 1:(nYear-1)]
  sSA <- S[2, 1:(nYear-1)]
  sAD <- matrix(0, nrow = nAge, ncol = nYear-1)

  if(nAgeC.S == 6){
    for(t in 1:(nYear-1)){
      sAD[2, t] <- S[3, t]
      for(a in 3:6) sAD[a, t] <- S[4, t] # prime-aged
      for(a in 7:9) sAD[a, t] <- S[5, t] # pre-senescent
      for(a in 10:nAge) sAD[a, t] <- S[6, t] # senescent
    }
  }else if(nAgeC.S == 12){
    for(t in 1:(nYear-1)){
      for(a in 2:11) sAD[a, t] <- S[a+1, t] # other adults
      for(a in 12:nAge) sAD[a, t] <- S[13, t] # greybeards
    }
  }else if(nAgeC.S == 20){
    for(t in 1:(nYear-1)){
      for(a in 2:nAge) sAD[a, t] <- S[a+1, t] # adults
    }
  }
  
  # priors for age-specific birth rate
  BR <- matrix(0, nrow = nAge, ncol = nYear-1)

  if(nAgeC.R == 6){
    for(t in 1:(nYear-1)){
      for(a in 2:4) BR[a, t] <- Ba[a-1, t]
      for(a in 5:6) BR[a, t] <- Ba[4, t]
      for(a in 7:10) BR[a, t] <- Ba[5, t]
      for(a in 11:nAge) BR[a, t] <- Ba[6, t]
    }
  }else if(nAgeC.R == 12){
    for(t in 1:(nYear-1)){
      for(a in 2:11) BR[a, t] <- Ba[a-1, t]
      for(a in 12:nAge) BR[a, t] <- Ba[11, t]
    }
  }else if(nAgeC.R == 20){
    for(t in 1:(nYear-1)){
      for(a in 2:nAge) BR[a, t] <- Ba[a-1, t]
    }
  }

  # priors for reproductive success
  sPY <- matrix(0, nrow = nAge, ncol = nYear-1)
  
  if(nAgeC.R == 6){
    for(t in 1:(nYear-1)){
      for(a in 2:4) sPY[a, t] <- Ra[a-1, t]
      for(a in 5:6) sPY[a, t] <- Ra[4, t]
      for(a in 7:10) sPY[a, t] <- Ra[5, t]
      for(a in 11:nAge) sPY[a, t] <- Ra[6, t]
    }
  }else if(nAgeC.R == 12){
    for(t in 1:(nYear-1)){
      for(a in 2:11) sPY[a, t] <- Ra[a-1, t]
      for(a in 12:nAge) sPY[a, t] <- Ra[11, t]
    }
  }else if(nAgeC.R == 20){
    for(t in 1:(nYear-1)){
      for(a in 2:nAge) sPY[a, t] <- Ra[a-1, t]
    }
  }
  
  ## Simulate observation parameters -------------------------------------------
  
  O <- runif(nYear, 0.1, 0.9)
  Mu.O <- runif(1, 0.1, 0.9)
  EpsilonT.O <- rnorm(nYear, 0, 0.2)
  SigmaT.O <- runif(1, 0.01, 2)
  
  
  ## Simulate initial population sizes -----------------------------------------
  
  # Actual numbers in 2008:
  # 5 female YFs in Sept, 6 SA1s, 5 SA2s, 21 adults
  # Wendy estimated 22.6% of the population was marked
  nYF        <- c(5*5, rep(NA, times = nYear-1))
  nYFa       <- matrix(1, nrow = nAge, ncol = nYear)
  nYFa[1:2,] <- 0 # too young to already have a YAF
  nSA        <- c(6*5, rep(NA, times = nYear-1))
  nAD        <- matrix(NA, nrow = nAge, ncol = nYear)
  nAD[,1]    <- c(0, 5*5, rep(2*5, times = 8), rep(1*5, times = nAge-10))
  nTOT       <- c(nYF[1] + nSA[1] + sum(nAD[2:nAge, 1]), rep(NA, times = nYear-1))
  
  # priors for initial population sizes
  initN.YF <- nYF[1]
  initN.SA <- nSA[1]
  initN.AD <- nAD[,1]
  
  log.initN.YF <- log(initN.YF)
  log.initN.SA <- log(initN.SA)
  log.initN.AD <- log(initN.AD)
  
  for(t in 1:(nYear-1)){
    # survival & birthdays
    nSA[t+1] <- pmax(10, rbinom(1, nYF[t], sYF[t]))
    nAD[2, t+1] <- pmax(10, rbinom(1, nSA[t], sSA[t]))
    for(a in 3:nAge) nAD[a, t+1] <- pmax(5, rbinom(1, nAD[a-1, t], sAD[a-1, t]))
    
    # then reproductive success
    for(a in 3:nAge) nYFa[a, t+1] <- pmax(1, rbinom(1, nAD[a-1, t], 0.5 * BR[a-1, t] * sPY[a-1, t]))
    nYF[t+1] <- sum(nYFa[3:nAge, t+1]) # total number of female YAFs every year
    nTOT[t+1] <- nYF[t+1] + nSA[t+1] + sum(nAD[2:nAge, t+1])
  }
  
  area <- rep(76.2, nYear)
  
  
  ## Assemble myinits list -----------------------------------------------------
  
  initList <- list(
    dens = dens,
    veg = veg,
    propF = propF,
    state = state,
    
    dens.true = dens.true,
    dens.cov = dens.cov,
    veg.true = veg.true,
    
    XiI.R = XiI.R,
    XiT.R = XiT.R,
    XiT.B = XiT.B,
    SigmaI.R = SigmaI.R,
    SigmaT.R = SigmaT.R,
    SigmaT.B = SigmaT.B,
    EpsilonI.R = EpsilonI.R,
    EpsilonT.R = EpsilonT.R,
    EpsilonT.B = EpsilonT.B,
    
    Mu.B = Mu.B,
    Mu.R = Mu.R,
    Bi = Bi,
    Ba = Ba,
    Ri = Ri,
    Ra = Ra,
    
    Mu.S = Mu.S,
    S = S,
    BR = BR,
    sPY = sPY,
    sYF = sYF,
    sSA = sSA,
    sAD = sAD,
    
    O = O,
    Mu.O = Mu.O,
    SigmaT.O = SigmaT.O,
    EpsilonT.O = EpsilonT.O,
    
    nYF = nYF,
    nYFa = nYFa,
    nSA = nSA,
    nAD = nAD,
    nTOT = nTOT,
    area = area,
    
    initN.YF = initN.YF,
    initN.SA = initN.SA,
    initN.AD = initN.AD,
    log.initN.YF = log.initN.YF,
    log.initN.SA = log.initN.SA,
    log.initN.AD = log.initN.AD
  )
  
  # append covariates
  if(envEffectsS){
    if(splitCovs == 3){
      initList <- c(initList, list(
        BetaD.Sy = BetaD.Sy,
        BetaD.Sp = BetaD.Sp,
        BetaD.So = BetaD.So,
        BetaV.Sy = BetaV.Sy,
        BetaV.Sp = BetaV.Sp,
        BetaV.So = BetaV.So))
    }else if(splitCovs == 2){
      initList <- c(initList, list(
        BetaD.Sy = BetaD.Sy,
        BetaD.So = BetaD.So,
        BetaV.Sy = BetaV.Sy,
        BetaV.So = BetaV.So))
    }else if(splitCovs == 1){
      initList <- c(initList, list(
        BetaD.S = BetaD.S,
        BetaV.S = BetaV.S))
    }
  }
  
  if(envEffectsR){
    initList <- c(initList, list(
      BetaD.R = BetaD.R))
  }
  
  # append random effects
  if(splitREs == 3){
    initList <- c(initList, list(
      XiT.Sy = XiT.Sy,
      XiT.Sp = XiT.Sp,
      XiT.So = XiT.So,
      SigmaT.Sy = SigmaT.Sy,
      SigmaT.Sp = SigmaT.Sp,
      SigmaT.So = SigmaT.So,
      EpsilonT.Sy = EpsilonT.Sy,
      EpsilonT.Sp = EpsilonT.Sp,
      EpsilonT.So = EpsilonT.So))
  }else if(splitREs == 2){
    initList <- c(initList, list(
      XiT.Sy = XiT.Sy,
      XiT.So = XiT.So,
      SigmaT.Sy = SigmaT.Sy,
      SigmaT.So = SigmaT.So,
      EpsilonT.Sy = EpsilonT.Sy,
      EpsilonT.So = EpsilonT.So))
  }else if(splitREs == 1){
    initList <- c(initList, list(
      XiT.S = XiT.S,
      SigmaT.S = SigmaT.S,
      EpsilonT.S = EpsilonT.S))
  }
  
  return(initList)
  
}

