#' Simulate initial values for IPM
#'
#' @param nR integer. Number of events in the reproductive success model. nR = 0 by default.
#' @param nID.S integer. Number of unique kangaroos in the survival model. nID.S = 0 by default.
#' @param nID.R integer. Number of unique kangaroos in the reproductive success model. nID.R = 0 by default.
#' @param nYear integer. Number of time steps in the model. nYear = 17 by default.
#' @param nAge integer. Number of ages, or maximum age, in the model. nAge = 17 by default.
#' @param nAgeC integer. Number of age classes in the model. nAgeC = 5 by default.
#' @param age.R vector of length nR of age of individuals in the analysis.
#' @param dens vector of lenth nYear of population density data.
#' @param veg vector of lenth nYear of available vegetation data.
#' @param win vector of lenth nYear of winter severity data.
#' @param nNoAge integer. Number of individuals of unknown age in the analysis. nNoAge = 0 by default.
#' @param nNoDens integer. Number of years when population density is unknown. nNoDens = 0 by default.
#' @param nNoVeg integer. Number of years when available vegetation is unknown. nNoVeg = 0 by default.
#' @param nNoWin integer. Number of years when winter severity is unknown. nNoWin = 0 by default.
#'
#' @returns list containing all initial values needed for the IPM.
#' @export
#'
#' @examples

simulateInits <- function(nR = 0, nID.S = 0, nID.R = 0, nYear = 17, nAge = 19, nAgeC = 5,
                          age.R, dens, veg, win, nNoAge = 0, nNoDens = 0, nNoVeg = 0, nNoWin = 0){
  
  # # for testing purposes
  # library(tidyverse)
  # library(readxl)
  # 
  # source("wrangleData_en.R")
  # enData <- wrangleData_en(dens.data = "data/WPNP_Methods_Results_January2025.xlsx",
  #                          veg.data  = "data/biomass data April 2009 - Jan 2025_updated Feb2025.xlsx",
  #                          wea.data  = "data/Prom_Weather_2008-2023_updated Jan2025 RB.xlsx",
  #                          wind.data = "data/POWER_Point_Daily_20080101_20241231_10M.csv")
  # 
  # source("wrangleData_rs.R")
  # rsData <- wrangleData_rs(rs.data = "data/RSmainRB_Mar25.xlsx",
  #                          obs.data = "data/PromObs_2008-2019.xlsx",
  #                          known.age = TRUE, cum.surv = TRUE, surv.sep1 = TRUE)
  # 
  # source("wrangleData_sv.R")
  # svData <- wrangleData_sv(surv.data = "data/PromSurvivalOct24.xlsx",
  #                          yafs.data = "data/RSmainRB_Mar25.xlsx")
  # 
  # nR <- rsData$nR
  # nID.S <- svData$nID
  # nID.R <- rsData$nID
  # nYear <- 17
  # nAge <- 19
  # nAgeC <- 5
  # age.R <- rsData$age.R
  # dens <- enData$dens
  # veg <- enData$veg
  # win <- enData$win
  # nNoAge <- svData$nNoAge
  # nNoDens <- enData$nNoDens
  # nNoVeg <- enData$nNoVeg
  # nNoWin <- enData$nNoWin
  
  if(missing(age.R)){
    age.R <- integer(nR)
  }
  
  
  ## Simulate latent states for input data -------------------------------------
  
  ## Survival model
  # true environment
  # dens.hat <- ifelse(is.na(dens), rnorm(length(dens), 0, .1), dens)
  # veg.hat  <- ifelse(is.na(veg), rnorm(length(veg), 0, .1), veg)
  
  dens.hat <- rnorm(nYear-1, 0, 1)
  veg.hat <- rnorm(nYear-1, 0, 1)
  
  # unobserved ages
  ageM <- sample(3:8, size = nNoAge, replace = T)
  
  # latent states
  Mu.Sp <- matrix(runif(nID.S * nYear, 0, 1), nrow = nID.S, ncol = nYear)
  Mu.Op <- matrix(runif(nID.S * nYear, 0, 1), nrow = nID.S, ncol = nYear)
  
  
  ## Simulate vital rate covariate effects -------------------------------------
  
  ## Survival model
  BetaA.S <- rnorm(nAgeC, 0, 1)
  BetaD.S <- rnorm(nAgeC, 0, 1)
  BetaV.S <- rnorm(nAgeC, 0, 1)
  # BetaDV.S <- rnorm(nAgeC, 0, 1)
  
  ## Reproductive success model
  # BetaD.rs <- rnorm(1, 0, 1)
  # BetaV.rs <- rnorm(1, 0, 1)
  # BetaW.rs <- rnorm(1, 0, 1)
  
  
  ## Simulate vital rate random effects ----------------------------------------
  
  ## Survival model
  # variance-covariance matrix
  zero  <- rep(0, nAgeC)
  Xi.S <- rnorm(nAgeC, 1, 0.1)
  
  Epsilon.S <- matrix(rnorm((nYear-1)*nAgeC, 0, 0.1),
                       ncol = nAgeC, nrow = (nYear-1))
  
  Gamma.S <- matrix(NA, ncol = nAgeC, nrow = nYear-1)
  
  for(t in 1:(nYear-1)){
    for(a in 1:nAgeC){
      Gamma.S[t,a] <- Xi.S[a] * Epsilon.S[t,a]
    }
  }
  
  Tau.S <- diag(nAgeC) + rnorm(nAgeC^2, 0, 0.1)
  Tau.S <- matrix(Tau.S, nrow = nAgeC)
  Tau.S <- (Tau.S + t(Tau.S)) / 2
  
  ## Reproductive success model
  EpsilonI.Ri <- rnorm(nID.R, 0, 1)
  EpsilonT.Ri <- rnorm(nYear, 0, 1)
  EpsilonT.Ra <- rnorm(nYear, 0, 1)
  
  SigmaI.Ri <- runif(1, 0, 10)
  SigmaT.Ri <- runif(1, 0, 10)
  SigmaT.Ra <- runif(1, 0, 10)

  
  ## Simulate yearly vital rates -----------------------------------------------
  
  ## Survival model
  # age-class-specific survival
  S <- matrix(NA, nrow = nAgeC, ncol = nYear-1)
  
  for(a in 1:nAgeC){
    for(t in 1:(nYear-1)){
      S[a, t] <- plogis(
        BetaA.S[a] +
        # BetaD.S[a] * dens.hat[t] +
        # BetaV.S[a] * veg.hat[t] +
        # BetaDV.S[a] * (dens.hat[t] * veg.hat[t]) +
        # BetaVR.S[a] * (veg.hat[t] / dens.hat[t]) +
        Gamma.S[t, a])
    }
  }

  ## Reproductive success model
  # age-specific reproductive success
  Mu.Ri <- runif(nAge, 0, 1)
  Mu.Ra <- runif(nAge, 0, 1)
  
  Ri <- numeric(nR)
  for(x in 1:nR) {
    Ri[x] <- plogis(
      qlogis(Mu.Ri[age.R[x]]) +
        # BetaD.rs * dens[year[x]] +
        # BetaV.rs * veg[year[x]] +
        # BetaW.rs * win[year[x]] +
        EpsilonI.Ri[x] +
        EpsilonT.Ri[x]
    )
  }
  
  Ra <- matrix(0, nrow = nAge, ncol = nYear)
  for(a in 1:nAge) {
    for(t in 1:nYear) {
      Ra[a, t] <- plogis(
        qlogis(Mu.Ra[a]) +
          # BetaD.rs * dens[t] +
          # BetaV.rs * veg[t] +
          # BetaW.rs * win[t] +
          EpsilonT.Ra[t]
      )
    }
  }
  
  ## Population model
  # breeding rate
  B <- runif(nYear, 0.5, 1) # raw means span 0.58-0.92
  
  # survival of PYs
  # to 1st Sept 1 when they become YAFs
  # sPY <- runif(nYear-1, 0.1, 1) # raw means span 0.24-0.95
  
  # Survival of YAFs
  # to 2nd Sept 1 when they become SA1s
  sYAF <- S[1, 1:(nYear-1)] # raw means span 0.01-0.88
  
  # Survival of SA1s to SA2 & SA2 to AD3
  # 2nd age class in our published CJS model
  sSA <- rbind(S[2, 1:(nYear-1)], S[2, 1:(nYear-1)])
  
  # Survival of all ADs
  sAD <- matrix(0, nrow = nAge, ncol = nYear-1)
  
  for(a in 3:6){
    sAD[a, 1:(nYear-1)] <- S[3, 1:(nYear-1)]
  }
  
  for(a in 7:9){
    sAD[a, 1:(nYear-1)] <- S[4, 1:(nYear-1)]
  }
  
  for(a in 10:nAge){
    sAD[a, 1:(nYear-1)] <- S[5, 1:(nYear-1)]
  }
  
  
  ## Simulate observation parameters -------------------------------------------
  
  # Recapture probabilities (CJS)
  O <- runif(nYear, 0.1, 0.9)
  Mu.O <- runif(1, 0.1, 0.9)
  Epsilon.O <- rnorm(nYear, 0, 0.2)
  Sigma.O <- runif(1, 0.01, 2) # or rnorm(1, 0.2, 0.1)
  
  
  ## Simulate initial population sizes -----------------------------------------
  
  # Actual numbers in 2008:
  # 5 female YAFs in Sept, 6 SA1s, 5 SA2s, 21 adults
  # Wendy estimated 22.6% of the population was marked
  
  nYAF    <- c(5*5, rep(NA, times = nYear-1))
  nYAFa   <- matrix(0, nrow = nAge, ncol = nYear)
  
  nSA     <- matrix(NA, nrow = 2, ncol = nYear)
  nSA[,1] <- c(6*5, 5*5)
  
  nAD     <- matrix(0, nrow = nAge, ncol = nYear)
  nAD[,1] <- c(0, 0, rep(2*5, times = 8), rep(1*5, times = nAge-10))
  
  nTOT   <- c(nYAF[1] + sum(nSA[1:2,1]) + sum(nAD[3:nAge,1]),
              rep(NA, times = nYear-1))
  
  for(t in 1:(nYear-1)){
    # survival & birthdays
    nSA[1, t+1] <- rbinom(1, nYAF[t], sYAF[t])
    nSA[2, t+1] <- rbinom(1, nSA[1, t], sSA[1, t])
    nAD[3, t+1] <- rbinom(1, nSA[2, t], sSA[2, t])
    for(a in 4:nAge){
      nAD[a, t+1] <- rbinom(1, nAD[a-1, t], sAD[a-1, t])
    }
    # then reproduction
    for(a in 3:nAge){
      nYAFa[a, t+1] <- rbinom(1, nAD[a, t+1], B[t+1] * Ra[a, t+1])
    }
    nYAF[t+1] <- sum(nYAFa[3:nAge, t+1])
    nTOT[t+1] <- nYAF[t+1] + sum(nSA[1:2, t+1]) + sum(nAD[3:nAge, t+1])
  }
  
  
  ## Assemble myinits list -----------------------------------------------------
  
  return(list(dens.hat = dens.hat,
              veg.hat = veg.hat,
              ageM = ageM,
              
              Mu.Sp = Mu.Sp,
              Mu.Op = Mu.Op,
              
              BetaA.S = BetaA.S,
              BetaD.S = BetaD.S,
              BetaV.S = BetaV.S,
              
              Xi.S = Xi.S,
              Epsilon.S = Epsilon.S,
              Gamma.S = Gamma.S,
              Tau.S = Tau.S,
              
              O = O,
              Mu.O = Mu.O,
              Epsilon.O = Epsilon.O,
              Sigma.O = Sigma.O,
              
              EpsilonI.Ri = EpsilonI.Ri,
              EpsilonT.Ri = EpsilonT.Ri,
              EpsilonT.Ra = EpsilonT.Ra,
              
              SigmaI.Ri = SigmaI.Ri,
              SigmaT.Ri = SigmaT.Ri,
              SigmaT.Ra = SigmaT.Ra,
              
              S = S,
              Mu.Ri = Mu.Ri,
              Mu.Ra = Mu.Ra,
              
              B = B,
              Ra = Ra,
              sYAF = sYAF,
              sSA = sSA,
              sAD = sAD,
              
              nYAF = nYAF,
              nSA = nSA,
              nAD = nAD,
              nTOT = nTOT
              ))
  
}

