#' Simulate initial values for IPM
#'
#' @param n integer. Number of events in the reproductive success model. n = 0 by default.
#' @param nID.S integer. Number of unique kangaroos in the survival model. nID.S = 0 by default.
#' @param nID.R integer. Number of unique kangaroos in the reproductive success model. nID.R = 0 by default.
#' @param nYear integer. Number of time steps in the model. nYear = 17 by default.
#' @param nAge integer. Number of ages, or maximum age, in the model. nAge = 17 by default.
#' @param nAgeC integer. Number of age classes in the model. nAgeC = 5 by default.
#' @param dens vector of yearly population density data.
#' @param veg vector of yearly available vegetation data.
#' @param win vector of yearly winter severity data.
#' @param nNoAge integer. Number of individuals of unknown age in the analysis. nNoAge = 0 by default.
#' @param nNoDens integer. Number of years when population density is unknown. nNoDens = 0 by default.
#' @param nNoVeg integer. Number of years when available vegetation is unknown. nNoVeg = 0 by default.
#' @param nNoWin integer. Number of years when winter severity is unknown. nNoWin = 0 by default.
#'
#' @returns list containing all initial values needed for the IPM.
#' @export
#'
#' @examples

simulateInits <- function(n = 0, nID.S = 0, nID.R = 0, nYear = 17, nAge = 17, nAgeC = 5,
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
  # n <- rsData$n
  # nID.S <- svData$nID
  # nID.R <- rsData$nID
  # nYear <- 17
  # nAge <- 17
  # nAgeC <- 5
  # dens <- enData$dens
  # veg <- enData$veg
  # win <- enData$win
  # nNoAge <- svData$nNoAge
  # nNoDens <- enData$N.noDens
  # nNoVeg <- enData$nNoVeg
  # nNoWin <- enData$nNoWin
  
  
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
  Mu.sp <- matrix(runif(nID.S * nYear, 0, 1), nrow = nID.S, ncol = nYear)
  Mu.op <- matrix(runif(nID.S * nYear, 0, 1), nrow = nID.S, ncol = nYear)
  
  
  ## Simulate vital rate covariate effects -------------------------------------
  
  ## Survival model
  BetaA.sv <- rnorm(nAgeC, 0, 1)
  BetaD.sv <- rnorm(nAgeC, 0, 1)
  BetaV.sv <- rnorm(nAgeC, 0, 1)
  # BetaDV.sv <- rnorm(nAgeC, 0, 1)
  
  ## Reproductive success model
  # BetaD.rs <- rnorm(1, 0, 1)
  # BetaV.rs <- rnorm(1, 0, 1)
  # BetaW.rs <- rnorm(1, 0, 1)
  
  
  ## Simulate vital rate random effects ----------------------------------------
  
  ## Survival model
  # variance-covariance matrix
  zero  <- rep(0, nAgeC)
  Xi.sv <- rnorm(nAgeC, 1, 0.1)
  
  Epsilon.sv <- matrix(rnorm((nYear-1)*nAgeC, 0, 0.1),
                       ncol = nAgeC, nrow = (nYear-1))
  
  Gamma.sv <- matrix(NA, ncol = nAgeC, nrow = nYear-1)
  
  for (t in 1:(nYear-1)){
    for (a in 1:nAgeC){
      Gamma.sv[t,a] <- Xi.sv[a] * Epsilon.sv[t,a]
    }
  }
  
  Tau.sv <- diag(nAgeC) + rnorm(nAgeC^2, 0, 0.1)
  Tau.sv <- matrix(Tau.sv, nrow = nAgeC)
  Tau.sv <- (Tau.sv + t(Tau.sv)) / 2
  
  ## Reproductive success model
  EpsilonI.rsI <- rnorm(nID.R, 0, 1)
  EpsilonT.rsI <- rnorm(nYear, 0, 1)
  EpsilonT.rsA <- rnorm(nYear, 0, 1)
  
  SigmaI.rsI <- runif(1, 0, 10)
  SigmaT.rsI <- runif(1, 0, 10)
  SigmaT.rsA <- runif(1, 0, 10)

  
  ## Simulate yearly vital rates -----------------------------------------------
  
  ## Survival model
  # age-class-specific survival
  sv <- matrix(NA, nrow = nAgeC, ncol = nYear-1)
  
  for(a in 1:nAgeC){
    for(t in 1:(nYear-1)){
      sv[a, t] <- plogis(
        BetaA.sv[a] +
        # BetaD.sv[a] * dens.hat[t] +
        # BetaV.sv[a] * veg.hat[t] +
        # BetaDV.sv[a] * (dens.hat[t] * veg.hat[t]) +
        # BetaVR.sv[a] * (veg.hat[t] / dens.hat[t]) +
        Gamma.sv[t, a])
    }
  }

  ## Reproductive success model
  # age-specific reproductive success
  Mu.rsI <- runif(nAge, 0, 1)
  Mu.rsA <- runif(nAge, 0, 1)
  
  rsI <- numeric(n)
  for (x in 1:n) {
    rsI[x] <- plogis(
      qlogis(Mu.rsI[age.R[x]]) +
        # BetaD.rs * dens[year[x]] +
        # BetaV.rs * veg[year[x]] +
        # BetaW.rs * win[year[x]] +
        EpsilonI.rsI[x] +
        EpsilonT.rsI[x]
    )
  }
  
  rsA <- numeric(nAge * nYear)
  for (a in 1:nAge) {
    for (t in 1:nYear) {
      idx <- (a - 1) * nYear + t
      rsA[idx] <- plogis(
        qlogis(Mu.rsA[a]) +
          # BetaD.rs * dens[t] +
          # BetaV.rs * veg[t] +
          # BetaW.rs * win[t] +
          EpsilonT.rsA[t]
        )
    }
  }
  
  ## Population model
  # breeding rate
  b <- runif(nYear-1, 0.5, 1) # raw means span 0.58-0.92
  
  # survival of PYs
  # to 1st Sept 1 when they become YAFs
  # svPY <- runif(nYear-1, 0.1, 1) # raw means span 0.24-0.95
  
  # TODO: DISCUSS NEW svPY HERE AS WELL!!
  for(a in 3:(nAge+2)){
    svPY[a, 1:(nYear-1)] <- rsA[a-2, 1:(nYear-1)] # age.R in RS model is age-2!
  }
  
  # Survival of YAFs
  # to 2nd Sept 1 when they become SA1s
  svYAF <- sv[1, 1:(nYear-1)] # raw means span 0.01-0.88
  
  # Survival of SA1s to SA2 & SA2 to AD3
  # 2nd age class in our published CJS model
  svSA <- rbind(sv[2, 1:(nYear-1)], sv[2, 1:(nYear-1)])
  
  # Survival of all ADs
  svAD <- matrix(0, nrow = nAge+2, ncol = nYear-1)
  
  for(a in 3:6){
    svAD[a, 1:(nYear-1)] <- sv[3, 1:(nYear-1)]
  }
  
  for(a in 7:9){
    svAD[a, 1:(nYear-1)] <- sv[4, 1:(nYear-1)]
  }
  
  for(a in 10:(nAge+2)){
    svAD[a, 1:(nYear-1)] <- sv[5, 1:(nYear-1)]
  }
  
  
  ## Simulate observation parameters -------------------------------------------
  
  # Recapture probabilities (CJS)
  ob <- runif(nYear, 0.1, 0.9)
  Mu.ob <- runif(1, 0.1, 0.9)
  Epsilon.ob <- rnorm(nYear, 0, 0.2)
  Sigma.ob <- runif(1, 0.01, 2) # or rnorm(1, 0.2, 0.1)
  
  
  ## Simulate initial population sizes -----------------------------------------
  
  # Actual numbers in 2008:
  # 5 female YAFs in Sept, 6 SA1s, 5 SA2s, 21 adults
  # Wendy estimated 22.6% of the population was marked
  
  nYAF    <- c(5*5, rep(NA, times = nYear-1)); nYAF
  nSA     <- matrix(NA, nrow = 2, ncol = nYear)
  nSA[,1] <- c(6*5, 5*5); nSA
  
  nAD     <- matrix(0, nrow = nAge+2, ncol = nYear); nAD
  nAD[,1] <- c(0, 0, rep(2*5, times = 8), rep(1*5, times = nAge-8)); nAD
  
  nTOT   <- c(nYAF[1] + sum(nSA[1:2,1]) + sum(nAD[3:(nAge+2),1]),
              rep(NA, times = nYear-1)); nTOT
  
  # TODO: DISCUSS NEW svPY HERE AS WELL!!
  for (t in 1:(nYear-1)){
    # nYAF[t+1]   <- rbinom(1, sum(nAD[3:(nAge+2), t]), b[t] * svPY[t])
    nYAF[t+1]   <- rbinom(1,
                          sum(nAD[3:(nAge+2), t]), 
                          sum(b[t] * svPY[3:(nAge+2), t] * nAD[3:(nAge+2), t]) / sum(nAD[3:(nAge+2), t]))
    
    nSA[1, t+1] <- rbinom(1, nYAF[t], svYAF[t])
    nSA[2, t+1] <- rbinom(1, nSA[1, t], svSA[1, t])
    nAD[3, t+1] <- rbinom(1, nSA[2, t], svSA[2, t])
    for (a in 4:(nAge+2)){
      nAD[a, t+1] <- rbinom(1, nAD[a-1, t], svAD[a-1, t])
    }
    nTOT[t+1] <- nYAF[t+1] + sum(nSA[1:2, t+1]) + sum(nAD[3:(nAge+2), t+1])
  }
  
  
  ## Assemble myinits list -----------------------------------------------------
  
  # TODO: UPDATE
  return(list(dens.hat = dens.hat,
              veg.hat = veg.hat,
              ageM = ageM,
              
              Mu.sp = Mu.sp,
              Mu.op = Mu.op,
              
              BetaA.sv = BetaA.sv,
              BetaD.sv = BetaD.sv,
              BetaV.sv = BetaV.sv,
              
              Xi.sv = Xi.sv,
              Epsilon.sv = Epsilon.sv,
              Gamma.sv = Gamma.sv,
              Tau.sv = Tau.sv,
              
              ob = ob,
              Mu.ob = Mu.ob,
              Epsilon.ob = Epsilon.ob,
              Sigma.ob = Sigma.ob,
              
              EpsilonI.rsI = EpsilonI.rsI,
              EpsilonT.rsI = EpsilonT.rsI,
              EpsilonT.rsA = EpsilonT.rsA,
              
              SigmaI.rsI = SigmaI.rsI,
              SigmaT.rsI = SigmaT.rsI,
              SigmaT.rsA = SigmaT.rsA,
              
              sv = sv,
              Mu.rsI = Mu.rsI,
              Mu.rsA = Mu.rsA,
              
              b = b,
              svPY = svPY,
              svYAF = svYAF,
              svSA = svSA,
              svAD = svAD,
              
              nYAF = nYAF,
              nSA = nSA,
              nAD = nAD,
              nTOT = nTOT
              ))
  
}

