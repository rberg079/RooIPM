#' Simulate initial values for IPM
#'
#' @param nYear integer. Number of time steps in the model. nYear = 17 by default.
#' @param nAge integer. Number of ages, or maximum age, in the model. nAge = 19 by default.
#' @param nR integer. Number of events in the reproductive success model. nR = 0 by default.
#' @param nID.S integer. Number of unique kangaroos in the survival model. nID.S = 0 by default.
#' @param nID.R integer. Number of unique kangaroos in the reproductive success model. nID.R = 0 by default.
#' @param year.R vector of length nR of years in the reproductive success analysis.
#' @param id.R vector of length nR of IDs of individuals in the reproductive success analysis.
#' @param age.R vector of length nR of age of individuals in the reproductive success analysis.
#' @param ageC.S vector of age classes to assign to actual ages in survival model.
#' @param ageC.R vector of age classes to assign to actual ages in reproductive success model.
#' @param dens vector of length nYear of population density data.
#' @param veg vector of length nYear of available vegetation data.
#' @param win vector of length nYear of winter severity data.
#' @param propF vector of length nYear of proportion of observations belonging to females.
#' @param envEffectsR logical. If TRUE, environmental covariates are included in RS model.
#' @param envEffectsS logical. If TRUE, environmental covariates are included in CJS model.
#'
#' @returns a list containing all initial values needed for the IPM.
#' @export
#'
#' @examples

simulateInits <- function(nYear = 17, nAge = 19, nR = 0, nID.S = 0, nID.R = 0,
                          year.R = 0, id.R = 0, age.R = 0, ageC.S, ageC.R,
                          dens, veg, win, propF = 0,
                          envEffectsS = TRUE, envEffectsR = TRUE){
  
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
  # source("wrangleData_rs.R")
  # rsData <- wrangleData_rs(rs.data = "data/RSmainRB_Mar25.xlsx",
  #                          obs.data = "data/PromObs_2008-2023.xlsx",
  #                          known.age = TRUE, cum.surv = FALSE)
  # 
  # source("wrangleData_sv.R")
  # svData <- wrangleData_sv(surv.data = "data/PromSurvivalOct24.xlsx",
  #                          yafs.data = "data/RSmainRB_Mar25.xlsx")
  # 
  # ageC.S <- c(1,2,2,3,3,3,3,4,4,4, rep(5,30)); ageC.S
  # ageC.R <- c(seq(from = 1, to = 19, by = 1), rep(19, times = 21)); ageC.R
  # 
  # nYear <- rsData$nYear
  # nAge <- rsData$nAge
  # nR <- rsData$nR
  # nID.S <- svData$nID
  # nID.R <- rsData$nID
  # year.R <- rsData$year.R
  # id.R <- rsData$id.R
  # age.R <- rsData$age.R
  # dens <- enData$dens
  # veg <- enData$veg
  # win <- enData$win
  # propF <- enData$propF
  # envEffectsS <- TRUE
  # envEffectsR <- TRUE
  # 
  # if(missing(age.R))  age.R  <- integer(nR)
  # if(missing(year.R)) year.R <- integer(nR)
  # if(missing(id.R))   id.R   <- integer(nR)
  
  
  ## Simulate latent states for input data -------------------------------------

  ## Survival model
  # missing values
  dens <- ifelse(is.na(dens), rnorm(nYear-1, 0, .1), dens)
  veg <- ifelse(is.na(veg), rnorm(nYear-1, 0, .1), veg)
  win <- ifelse(is.na(win), rnorm(nYear-1, 0, .1), win)
  propF <- ifelse(is.na(propF), rnorm(nYear, .7, .05), propF)
  
  # true environment
  dens.hat <- ifelse(is.na(dens), rnorm(nYear-1, 0, .1), dens)
  veg.hat  <- ifelse(is.na(veg), rnorm(nYear-1, 0, .1), veg)
  win.hat  <- ifelse(is.na(win), rnorm(nYear-1, 0, .1), win)
  
  # latent states
  Mu.Sp <- matrix(runif(nID.S * nYear, 0, 1), nrow = nID.S, ncol = nYear)
  Mu.Op <- matrix(runif(nID.S * nYear, 0, 1), nrow = nID.S, ncol = nYear)
  
  
  ## Simulate vital rate covariate effects -------------------------------------
  
  nAgeC.S = max(ageC.S)
  nAgeC.R = max(ageC.R)
  
  ## Survival model
  BetaA.S <- c(rnorm(1, 1.0, 0.2),
               rnorm(1, 2.4, 0.2),
               rnorm(1, 2.8, 0.2),
               rnorm(1, 2.4, 0.2),
               rnorm(1, 1.0, 0.2))
  
  if(envEffectsS){
    BetaD.S <- runif(nAgeC.S, -5, 5)
    BetaV.S <- runif(nAgeC.S, -5, 5)
    BetaW.S <- runif(nAgeC.S, -5, 5)
  }
  
  ## Reproductive success model
  if(envEffectsR){
    BetaD.R <- runif(1, -2, 2)
    BetaV.R <- runif(1, -2, 2)
    BetaW.R <- runif(1, -2, 2)
  }
  
  
  ## Simulate vital rate random effects ----------------------------------------
  
  ## Survival model
  # variance-covariance matrix
  zero  <- rep(0, nAgeC.S)
  Xi.S <- rnorm(nAgeC.S, 1, 0.1)
  
  Epsilon.S <- matrix(rnorm((nYear-1)*nAgeC.S, 0, 0.1),
                       ncol = nAgeC.S, nrow = (nYear-1))
  
  Gamma.S <- matrix(NA, ncol = nAgeC.S, nrow = nYear-1)
  
  for(t in 1:(nYear-1)){
    for(a in 1:nAgeC.S){
      Gamma.S[t,a] <- Xi.S[a] * Epsilon.S[t,a]
    }
  }
  
  Tau.S <- diag(nAgeC.S) + rnorm(nAgeC.S^2, 0, 0.1)
  Tau.S <- matrix(Tau.S, nrow = nAgeC.S)
  Tau.S <- (Tau.S + t(Tau.S)) / 2
  
  ## Reproductive success model
  EpsilonI.R <- rnorm(nID.R, 0, 1)
  EpsilonT.R <- rnorm(nYear-1, 0, 1)
  EpsilonT.B <- rnorm(nYear-1, 0, 1)

  SigmaI.R <- runif(1, 0, 5)
  SigmaT.R <- runif(1, 0, 5)
  SigmaT.B <- runif(1, 0, 5)

  
  ## Simulate yearly vital rates -----------------------------------------------
  
  ## Survival model
  # age-class-specific survival
  S <- matrix(NA, nrow = nAgeC.S, ncol = nYear-1)
  
  for(a in 1:nAgeC.S){
    for(t in 1:(nYear - 1)){
      if(envEffectsS){
        S[a, t] <- plogis(
          BetaA.S[a] +
            BetaD.S[a] * dens.hat[t] +
            BetaV.S[a] * veg.hat[t] +
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
  Mu.R <- c(0, rep(runif(nAgeC.R-1, 0, 1)))

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
          BetaD.R * dens.hat[year.R[x]] +
          BetaV.R * veg.hat[year.R[x]] +
          BetaW.R * win.hat[year.R[x]] +
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
            BetaD.R * dens.hat[t] +
            BetaV.R * veg.hat[t] +
            BetaW.R * win.hat[t] +
            EpsilonT.R[t])
      }else{
        Ra[a, t] <- plogis(
          qlogis(Mu.R[a]) +
            EpsilonT.R[t])
      }
    }
  }
  
  ## Population model
  # survival of YAFs to 2nd Sept 1 when they become SA1s
  sYAF <- S[1, 1:(nYear-1)]
  
  # survival of SA1s to SA2 (now AD2)
  sSA <- S[2, 1:(nYear-1)]
  
  # survival of all ADs
  sAD <- matrix(0, nrow = nAge, ncol = nYear-1)
  sAD[2, 1:(nYear-1)] <- S[2, 1:(nYear-1)]
  
  for(a in 3:6){ # prime-aged
    sAD[a, 1:(nYear-1)] <- S[3, 1:(nYear-1)]
  }
  for(a in 7:9){ # pre-senescent
    sAD[a, 1:(nYear-1)] <- S[4, 1:(nYear-1)]
  }
  for(a in 10:nAge){ # senescent
    sAD[a, 1:(nYear-1)] <- S[5, 1:(nYear-1)]
  }
  
  # reproductive success
  rAD <- matrix(0, nrow = nAge, ncol = nYear-1)
  for(a in 1:nAgeC.R){
    rAD[a, 1:(nYear-1)] <- Ra[a, 1:(nYear-1)]
  }
  if(nAge > nAgeC.R){
    for(a in (nAgeC.R+1):nAge){
      rAD[a, 1:(nYear-1)] <- Ra[nAgeC.R, 1:(nYear-1)]
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
  # 5 female YAFs in Sept, 6 SA1s, 5 SA2s, 21 adults
  # Wendy estimated 22.6% of the population was marked
  
  nYAF  <- c(5*5, rep(NA, times = nYear-1)); nYAF
  nYAFa <- matrix(1, nrow = nAge, ncol = nYear)
  nYAFa[1:2,] <- 0; nYAFa
  
  nSA <- c(6*5, rep(NA, times = nYear-1)); nSA
  
  nAD     <- matrix(NA, nrow = nAge, ncol = nYear)
  nAD[,1] <- c(0, 5*5, rep(2*5, times = 8), rep(1*5, times = nAge-10)); nAD
  
  nTOT <- c(nYAF[1] + nSA[1] + sum(nAD[2:nAge, 1]), rep(NA, times = nYear-1)); nTOT
  
  for(t in 1:(nYear-1)){
    # survival & birthdays
    nSA[t+1] <- pmax(10, rbinom(1, nYAF[t], sYAF[t]))
    nAD[2, t+1] <- pmax(10, rbinom(1, nSA[t], sSA[t]))
    for(a in 3:nAge){
      nAD[a, t+1] <- pmax(5, rbinom(1, nAD[a-1, t], sAD[a-1, t]))
    }
    # then reproductive success
    for(a in 3:nAge){
      nYAFa[a, t+1] <- pmax(1, rbinom(1, nAD[a-1, t], 0.5 * Bt[t] * rAD[a-1, t]))
    }
    nYAF[t+1] <- sum(nYAFa[3:nAge, t+1])
    nTOT[t+1] <- nYAF[t+1] + nSA[t+1] + sum(nAD[2:nAge, t+1])
  }
  
  ab <- round(pmax((nTOT / pmax(propF, .4)) + rnorm(length(nTOT), 0, 2), 1))
  # pmax(propF, .4) returns .4 if propF falls below it
  # pmax(..., 1) returns 1 if ab falls below it
  # so propF is at least 40% & ab at least 1
  
  # nYAF; nSA; nAD; nTOT; ab
  
  ## Assemble myinits list -----------------------------------------------------
  
  initList <- list(
    dens = dens,
    veg = veg,
    win = win,
    propF = propF,
    dens.hat = dens.hat,
    veg.hat = veg.hat,
    win.hat = win.hat,
    
    Mu.Sp = Mu.Sp,
    Mu.Op = Mu.Op,
    
    Xi.S = Xi.S,
    Epsilon.S = Epsilon.S,
    Gamma.S = Gamma.S,
    Tau.S = Tau.S,
    
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
    rAD = rAD,
    
    S = S,
    sYAF = sYAF,
    sSA = sSA,
    sAD = sAD,
    
    O = O,
    Mu.O = Mu.O,
    EpsilonT.O = EpsilonT.O,
    SigmaT.O = SigmaT.O,
    
    nYAF = nYAF,
    nYAFa = nYAFa,
    nSA = nSA,
    nAD = nAD,
    nTOT = nTOT,
    ab = ab
  )
  
  if(envEffectsS){
    initList <- c(initList, list(
      BetaA.S = BetaA.S,
      BetaD.S = BetaD.S,
      BetaV.S = BetaV.S
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

