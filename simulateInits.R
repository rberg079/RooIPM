#' Simulate initial values for IPM
#'
#' @param nYear integer. Number of time steps in the model. nYear = 17 by default.
#' @param nAge integer. Number of ages, or maximum age, in the model. nAge = 19 by default.
#' @param ageClasses integer. Number of age classes to be considered. ageClasses = 20 by default.
#' @param nID.S integer. Number of unique kangaroos in the survival model. nID.S = 0 by default.
#' @param ageC.S integer vector. Age classes to assign to actual ages in survival model.
#' @param nR integer. Number of events in the reproductive success model. nR = 0 by default.
#' @param nID.R integer. Number of unique kangaroos in the reproductive success model. nID.R = 0 by default.
#' @param year.R integer vector. Year of each event in the reproductive success analysis.
#' @param id.R integer vector. Maternal ID of each event in the reproductive success analysis.
#' @param age.R integer vector. Maternal age of each event in the reproductive success analysis.
#' @param ageC.R integer vector. Age classes to assign to actual ages in reproductive success model.
#' @param dens vector of length nYear of population density data.
#' @param veg vector of length nYear of available vegetation data.
#' @param win vector of length nYear of winter severity data.
#' @param propF vector of length nYear of proportion of observations belonging to females.
#' @param envEffectsS logical. If TRUE, environmental covariates are included in CJS model.
#' @param envEffectsR logical. If TRUE, environmental covariates are included in RS model.
#'
#' @returns a list containing all initial values needed for the IPM.
#' @export
#'
#' @examples

simulateInits <- function(nYear = 17, nAge = 19, ageClasses = 20, nID.S = 0, ageC.S,
                          nR = 0, nID.R = 0, year.R = 0, id.R = 0, age.R = 0, ageC.R,
                          dens, veg, win, propF = 0, envEffectsS = TRUE, envEffectsR = TRUE){
  
  # # for testing purposes
  # library(readxl)
  # library(tidyverse)
  # 
  # source("wrangleData_en.R")
  # enData <- wrangleData_en(dens.data = "data/abundanceData_Proteus.csv",
  #                          veg.data  = "data/biomass data April 2009 - Jan 2025_updated Feb2025.xlsx",
  #                          wea.data  = "data/Prom_Weather_2008-2023_updated Jan2025 RB.xlsx",
  #                          wind.data = "data/POWER_Point_Daily_20080101_20241231_10M.csv",
  #                          obs.data  = "data/PromObs_2008-2023.xlsx",
  #                          list      = "data/PromlistAllOct24.xlsx")
  # 
  # source('wrangleData_sv.R')
  # svData <- wrangleData_sv(surv.data = "data/PromSurvivalOct24.xlsx",
  #                          yafs.data = "data/RSmainRB_Mar25.xlsx",
  #                          ageClasses = 6, known.age = TRUE)
  # 
  # source('wrangleData_rs.R')
  # rsData <- wrangleData_rs(rs.data = "data/RSmainRB_Mar25.xlsx",
  #                          obs.data = "data/PromObs_2008-2023.xlsx",
  #                          ageClasses = 6, known.age = TRUE, cum.surv = FALSE)
  # 
  # nYear <- rsData$nYear
  # nAge <- rsData$nAge
  # ageClasses <- 6
  # nID.S <- svData$nID.S
  # ageC.S <- svData$ageC.S
  # nAgeC.S <- svData$nAgeC.S
  # nR <- rsData$nR
  # nID.R <- rsData$nID.R
  # year.R <- rsData$year.R
  # id.R <- rsData$id.R
  # age.R <- rsData$age.R
  # ageC.R <- rsData$ageC.R
  # nAgeC.R <- rsData$nAgeC.R
  # dens <- enData$dens
  # veg <- enData$veg
  # win <- enData$win
  # propF <- enData$propF
  # envEffectsR <- TRUE
  # envEffectsS <- TRUE
  # 
  # if(missing(age.R))  age.R  <- integer(nR)
  # if(missing(year.R)) year.R <- integer(nR)
  # if(missing(id.R))   id.R   <- integer(nR)
  
  
  ## Simulate latent states for input data -------------------------------------
  
  ## Survival model
  # missing values
  nNoDens <- sum(is.na(dens))
  nNoVeg  <- sum(is.na(veg))
  nNoWin  <- sum(is.na(win))
  nNoProp <- sum(is.na(propF))
  
  dens <- round(ifelse(is.na(dens), rnorm(nYear, 3.9, .4), dens), 2)
  veg <- round(ifelse(is.na(veg), rnorm(nYear, 0, .1), veg), 4)
  win <- round(ifelse(is.na(win), rnorm(nYear, 0, .1), win), 4)
  propF <- round(ifelse(is.na(propF), pmax(pmin(rnorm(nYear, .7, .1), 0.99), 0.4), propF), 4)
  
  # true environment
  dens.true <- dens
  veg.true  <- veg
  win.true  <- win
  
  # latent states
  Mu.Sp <- matrix(runif(nID.S * nYear, 0, 1), nrow = nID.S, ncol = nYear)
  Mu.Op <- matrix(runif(nID.S * nYear, 0, 1), nrow = nID.S, ncol = nYear)
  
  
  ## Simulate vital rate covariate effects -------------------------------------
  
  nAgeC.S <- max(ageC.S)
  nAgeC.R <- max(ageC.R)
  
  ## Survival model
  if(ageClasses == 6){
    BetaA.S <- c(rnorm(1, 1.0, 0.2),
                 rnorm(1, 1.0, 0.2),
                 rnorm(1, 2.4, 0.2),
                 rnorm(1, 2.8, 0.2),
                 rnorm(1, 2.4, 0.2),
                 rnorm(1, 1.0, 0.2))
  }else if(ageClasses == 12){
    BetaA.S <- c(rnorm(2, 1.0, 0.2),
                 rnorm(1, 2.4, 0.2),
                 rnorm(4, 2.8, 0.2),
                 rnorm(3, 2.4, 0.2),
                 rnorm(3, 1.0, 0.2))
  }else if(ageClasses == 20){
    BetaA.S <- c(rnorm(1, 1.0, 0.2), rnorm(1, 1.0, 0.2),
                 rnorm(1, 2.4, 0.2), rnorm(1, 2.8, 0.2),
                 rnorm(1, 2.8, 0.2), rnorm(1, 2.8, 0.2),
                 rnorm(1, 2.8, 0.2), rnorm(1, 2.4, 0.2),
                 rnorm(1, 2.4, 0.2), rnorm(1, 2.4, 0.2),
                 rep(rnorm(1, 1.0, 0.2), 10)) 
  }
  
  if(envEffectsS){
    BetaD.S <- runif(nAgeC.S, -1, 1)
    BetaV.S <- runif(nAgeC.S, -1, 1)
    BetaW.S <- runif(nAgeC.S, -1, 1)
  }
  
  ## Reproductive success model
  if(envEffectsR){
    BetaD.R <- runif(1, -1, 1)
    BetaV.R <- runif(1, -1, 1)
    BetaW.R <- runif(1, -1, 1)
  }
  
  
  ## Simulate vital rate random effects ----------------------------------------
  
  ## Survival model
  # variance-covariance matrix
  Xi.S <- rnorm(nAgeC.S, 1, 0.1)
  
  Epsilon.S <- matrix(rnorm((nYear-1)*nAgeC.S, 0, 0.1),
                      nrow = (nYear-1), ncol = nAgeC.S)
  
  Gamma.S <- matrix(NA, ncol = nAgeC.S, nrow = nYear-1)
  
  for(t in 1:(nYear-1)){
    for(a in 1:nAgeC.S){
      Gamma.S[t,a] <- Xi.S[a] * Epsilon.S[t,a]
    }
  }
  
  # Tau.S <- diag(nAgeC.S) + rnorm(nAgeC.S^2, 0, 0.1)
  # Tau.S <- (Tau.S + t(Tau.S)) / 2
  
  # alternative init Tau.S that guarantees positive-definite values
  # (apparently potentially problematic the way I had it above)
  A <- matrix(rnorm(nAgeC.S^2, 0, 0.1), nAgeC.S, nAgeC.S)
  Tau.S <- crossprod(A) + diag(nAgeC.S)  # positive-definite
  
  ## Reproductive success model
  # latent standard normals
  XiI.R <- rnorm(nID.R, 0, 1)
  XiT.R <- rnorm(nYear-1, 0, 1)
  XiT.B <- rnorm(nYear-1, 0, 1)

  # scales
  SigmaI.R <- runif(1, .5, 2)
  SigmaT.R <- runif(1, .5, 2)
  SigmaT.B <- runif(1, .5, 2)
  
  # scaled random effects
  EpsilonI.R <- XiI.R * SigmaI.R
  EpsilonT.R <- XiT.R * SigmaT.R
  EpsilonT.B <- XiT.B * SigmaT.B

  
  ## Simulate yearly vital rates -----------------------------------------------
  
  ## Survival model
  # age-class-specific survival
  S <- matrix(NA, nrow = nAgeC.S, ncol = nYear-1)
  
  for(a in 1:nAgeC.S){
    for(t in 1:(nYear - 1)){
      if(envEffectsS){
        S[a, t] <- plogis(
          BetaA.S[a] +
            BetaD.S[a] * dens.true[t] +
            BetaV.S[a] * veg.true[t] +
            BetaW.S[a] * win.true[t] +
            Gamma.S[t, a])
      }else{
        S[a, t] <- plogis(
          BetaA.S[a] +
            Gamma.S[t, a])
      }
    }
  }

  # Reproductive success model
  # age-specific reproductive success
  Mu.B <- runif(1, 0.4, 1)
  Mu.R <- c(rep(runif(nAgeC.R, 0, 1)))

  Bi <- numeric(nR)
  for(x in 1:nR){
    Bi[x] <- plogis(
      qlogis(Mu.B) + EpsilonT.B[year.R[x]])
  }

  Bt <- numeric(nYear-1)
  for(t in 1:(nYear-1)){
    Bt[t] <- plogis(qlogis(Mu.B) + EpsilonT.B[t])
  }
  
  Ri <- numeric(nR)
  for(x in 1:nR){
    if(envEffectsR){
      Ri[x] <- plogis(
        qlogis(Mu.R[ageC.R[age.R[x]]]) +
          BetaD.R * dens.true[year.R[x]] +
          BetaV.R * veg.true[year.R[x]] +
          BetaW.R * win.true[year.R[x]] +
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
            BetaD.R * dens.true[t] +
            BetaV.R * veg.true[t] +
            BetaW.R * win.true[t] +
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
      for(a in 2:19) sAD[a, t] <- S[a+1, t] # adults
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
      for(a in 2:19) sPY[a, t] <- Ra[a-1, t]
    }
  }
  
  
  ## Simulate observation parameters -------------------------------------------
  
  # Recapture probabilities (CJS)
  O <- runif(nYear, 0.1, 0.9)
  Mu.O <- runif(1, 0.1, 0.9)
  EpsilonT.O <- rnorm(nYear, 0, 0.2)
  SigmaT.O <- runif(1, 0.01, 2) # or rnorm(1, 0.2, 0.1)
  
  
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
  
  for(t in 1:(nYear-1)){
    # survival & birthdays
    nSA[t+1] <- pmax(10, rbinom(1, nYF[t], sYF[t]))
    nAD[2, t+1] <- pmax(10, rbinom(1, nSA[t], sSA[t]))
    for(a in 3:nAge) nAD[a, t+1] <- pmax(5, rbinom(1, nAD[a-1, t], sAD[a-1, t]))
    
    # then reproductive success
    for(a in 3:nAge) nYFa[a, t+1] <- pmax(1, rbinom(1, nAD[a-1, t], 0.5 * Bt[t] * sPY[a-1, t]))
    nYF[t+1] <- sum(nYFa[3:nAge, t+1]) # total number of female YAFs every year
    nTOT[t+1] <- nYF[t+1] + nSA[t+1] + sum(nAD[2:nAge, t+1])
  }
  
  area <- rep(76.2, nYear)
  
  # ab <- round(pmax((nTOT / pmax(propF, .4)) + rnorm(length(nTOT), 0, 2), 1))
  # # pmax(propF, .4) returns .4 if propF falls below it
  # # pmax(..., 1) returns 1 if ab falls below it
  # # so propF is at least 40% & ab at least 1

  # nYF; nSA; nAD; nTOT; ab
  
  ## Assemble myinits list -----------------------------------------------------
  
  initList <- list(
    dens = dens,
    veg = veg,
    win = win,
    propF = propF,
    dens.true = dens.true,
    veg.true = veg.true,
    win.true = win.true,
    
    Mu.Sp = Mu.Sp,
    Mu.Op = Mu.Op,
    BetaA.S = BetaA.S,
    
    Xi.S = Xi.S,
    Epsilon.S = Epsilon.S,
    Gamma.S = Gamma.S,
    Tau.S = Tau.S,
    
    XiI.R = XiI.R,
    XiT.R = XiT.R,
    XiT.B = XiT.B,
    
    EpsilonI.R = EpsilonI.R,
    EpsilonT.R = EpsilonT.R,
    EpsilonT.B = EpsilonT.B,
    
    SigmaI.R = SigmaI.R,
    SigmaT.R = SigmaT.R,
    SigmaT.B = SigmaT.B,
    
    Mu.B = Mu.B,
    Mu.R = Mu.R,
    Bi = Bi,
    Bt = Bt,
    Ri = Ri,
    Ra = Ra,
    
    S = S,
    sPY = sPY,
    sYF = sYF,
    sSA = sSA,
    sAD = sAD,
    
    O = O,
    Mu.O = Mu.O,
    EpsilonT.O = EpsilonT.O,
    SigmaT.O = SigmaT.O,
    
    nYF = nYF,
    nYFa = nYFa,
    nSA = nSA,
    nAD = nAD,
    nTOT = nTOT,
    area = area
  )
  
  if(envEffectsS){
    initList <- c(initList, list(
      BetaD.S = BetaD.S,
      BetaV.S = BetaV.S,
      BetaW.S = BetaW.S
    ))
  }
  
  if(envEffectsR){
    initList <- c(initList, list(
      BetaD.R = BetaD.R,
      BetaV.R = BetaV.R,
      BetaW.R = BetaW.R
    ))
  }
  
  return(initList)
  
}

