#' Simulate initial values for IPM
#'
#' @param N
#' @param N.id 
#' @param N.year integer. Number of time steps in the model. N.year = 17 by default.
#' @param N.age integer. Number of adult age classes, or maximum age. N.age = 22 by default.
#' @param N.ageC integer. Number of age classes in CJS model. N.ageC = 5 by default.
#' @param dens vector. Where the density data is stored.
#' @param veg vector. Where the vegetation data is stored.
#' @param win 
#' @param N.noAge integer. Number of individuals of unknown age in the analysis.
#' @param N.noDens 
#' @param N.noVeg 
#' @param N.noWin 
#'
#' @returns list containing all initial values needed for the IPM.
#' @export
#'
#' @examples

simulateInits <- function(N, N.id, N.year = 17, N.age, N.ageC = 5, dens, veg, win,
                          N.noAge, N.noDens = 0, N.noVeg = 0, N.noWin = 0){
  
  # for testing purposes
  source("wrangleData_env.R")
  enData <- wrangleData_env(dens.data = "data/WPNP_Methods_Results_January2025.xlsx",
                            veg.data  = "data/biomass data April 2009 - Jan 2025_updated Feb2025.xlsx",
                            wea.data  = "data/Prom_Weather_2008-2023_updated Jan2025 RB.xlsx",
                            wind.data = "data/POWER_Point_Daily_20080101_20241231_10M.csv")
  
  source("wrangleData_rs.R")
  rsData <- wrangleData_rs(rs.data = "data/RSmainRB_Mar25.xlsx",
                           obs.data = "data/PromObs_2008-2019.xlsx",
                           known.age = TRUE, cum.surv = TRUE, surv.sep1 = TRUE)
  
  source("wrangleData_surv.R")
  svData <- wrangleData_surv(surv.data = "data/PromSurvivalOct24.xlsx",
                             yafs.data = "data/RSmainRB_Mar25.xlsx")
  
  N <- rs.dat$N
  N.id <- rs.dat$N.id
  N.age <- rs.dat$N.age
  dens <- en.dat$dens
  veg <- en.dat$veg
  win <- en.dat$win
  N.noAge <- surv$N.noAge
  N.noDens <- en.dat$N.noDens
  N.noVeg <- en.dat$N.noVeg
  N.noWin <- en.dat$N.noWin
  
  
  ## Simulate latent states for input data -------------------------------------
  
  ## Survival model
  # true environment
  dens.hat <- ifelse(is.na(dens), rnorm(length(dens), 0, .1), dens)
  veg.hat  <- ifelse(is.na(veg), rnorm(length(veg), 0, .1), veg)
  
  # unobserved ages
  ageM <- sample(3:8, size = N.noAge, replace = T)
  
  # latent states
  Mu.sp <- runif(1, 0, 1)
  Mu.op <- runif(1, 0, 1)
  
  
  ## Simulate vital rate covariate effects -------------------------------------
  
  ## Survival model
  BetaA.sv <- rnorm(N.ageC, 0, 1)
  BetaD.sv <- rnorm(N.ageC, 0, 1)
  BetaV.sv <- rnorm(N.ageC, 0, 1)
  # BetaDV.sv <- rnorm(N.ageC, 0, 1)
  
  ## Reproductive success model
  # BetaD.rs <- rnorm(1, 0, 1)
  # BetaV.rs <- rnorm(1, 0, 1)
  # BetaW.rs <- rnorm(1, 0, 1)
  
  
  ## Simulate vital rate random effects ----------------------------------------
  
  ## Survival model
  # variance-covariance matrix
  zero  <- rep(0, N.ageC)
  Xi.sv <- rnorm(N.ageC, 1, 0.1)
  
  Epsilon.sv <- matrix(rnorm((N.year-1)*N.ageC, 0, 0.1),
                    ncol = N.ageC, nrow = (N.year-1))
  
  Gamma.sv <- matrix(NA, ncol = N.ageC, nrow = N.year-1)
  
  for (t in 1:(N.year-1)){
    for (a in 1:N.ageC){
      Gamma.sv[t,a] <- Xi.sv[a] * Epsilon.sv[t,a]
    }
  }
  
  Tau.sv <- diag(N.ageC) + rnorm(N.ageC^2, 0, 0.1)
  Tau.sv <- matrix(Tau.sv, nrow = N.ageC)
  Tau.sv <- (Tau.sv + t(Tau.sv)) / 2
  
  ## Reproductive success model
  EpsilonI.rsI <- rnorm(N.id, 0, 1)
  EpsilonT.rsI <- rnorm(N.year, 0, 1)
  EpsilonT.rsA <- rnorm(N.year, 0, 1)
  
  SigmaI.rsI <- runif(1, 0, 10)
  SigmaT.rsI <- runif(1, 0, 10)
  SigmaT.rsA <- runif(1, 0, 10)

  
  ## Simulate yearly vital rates -----------------------------------------------
  
  ## Survival model
  # age-class-specific survival
  sv <- matrix(NA, nrow = N.ageC, ncol = N.year-1)
  
  for(a in 1:N.ageC){
    for(t in 1:(N.year-1)){
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
  # rsI <- numeric(N)
  # for (i in 1:N) {
  #   rsI[i] <- plogis(
  #     qlogis(Mu.rsI[age[i]]) +
  #       # BetaD.rs * dens[t] +
  #       # BetaV.rs * veg[t] +
  #       # BetaW.rs * win[t] +
  #       EpsilonI.rsI[id[i]] +
  #       EpsilonT.rsI[year[i]]
  #   )
  # }
  # 
  # rsA <- numeric(N.age * N.year)
  # for (a in 1:N.age) {
  #   for (t in 1:N.year) {
  #     idx <- (a - 1) * N.year + t
  #     rsA[idx] <- plogis(
  #       qlogis(Mu.rsA[a]) +
  #         # BetaD.rs * dens[t] +
  #         # BetaV.rs * veg[t] +
  #         # BetaW.rs * win[t] +
  #         EpsilonT.rsA[t]
  #       )
  #   }
  # }
  
  Mu.rsI <- runif(N.age, 0, 1)
  Mu.rsA <- runif(N.age, 0, 1)
  
  ## Population model
  # breeding rate
  b <- runif(N.year-1, 0.5, 1) # raw means span 0.58-0.92
  
  # survival of PYs
  # to 1st Sept 1 when they become YAFs
  s.PY <- runif(N.year-1, 0.1, 1) # raw means span 0.24-0.95
  
  # Survival of YAFs
  # to 2nd Sept 1 when they become SA1s
  s.YAF <- sv[1, 1:(N.year-1)] # raw means span 0.01-0.88
  
  # Survival of SA1s to SA2 & SA2 to AD3
  # 2nd age class in our published CJS model
  s.SA <- rbind(sv[2, 1:(N.year-1)], sv[2, 1:(N.year-1)])
  
  # Survival of all ADs
  s.AD <- matrix(0, nrow = N.age, ncol = N.year-1)
  
  for(a in 3:6){
    s.AD[a, 1:(N.year-1)] <- sv[3, 1:(N.year-1)]
  }
  
  for(a in 7:9){
    s.AD[a, 1:(N.year-1)] <- sv[4, 1:(N.year-1)]
  }
  
  for(a in 10:N.age){
    s.AD[a, 1:(N.year-1)] <- sv[5, 1:(N.year-1)]
  }
  
  
  ## Simulate observation parameters -------------------------------------------
  
  # Recapture probabilities (CJS)
  year.p <- rnorm(N.year, 0, 0.2)
  sd.p <- rnorm(1, 0.2, 0.1)
  
  
  ## Simulate initial population sizes -----------------------------------------
  
  # Actual numbers in 2008:
  # 5 female YAFs in Sept, 6 SA1s, 5 SA2s, 21 adults
  # Wendy estimated 22.6% of the population was marked
  
  YAF    <- c(5*5, rep(NA, times = N.year-1))
  SA     <- matrix(NA, nrow = 2, ncol = N.year)
  SA[,1] <- c(6*5, 5*5)
  
  AD     <- matrix(NA, nrow = N.age, ncol = N.year)
  AD[,1] <- c(0, 0, rep(2*5, times = (N.age-2)/2), rep(1*5, times = (N.age-2)/2))
  
  Ntot   <- c(YAF[1] + sum(SA[1:2,1]) + sum(AD[3:N.age,1]),
              rep(NA, times = N.year-1))
  
  
  ## Assemble myinits list -----------------------------------------------------
  
  # TODO: UPDATE
  return(list(dens.hat = dens.hat,
              veg.hat = veg.hat,
              ageM = ageM,
              
              BetaA.sv = BetaA.sv,
              BetaD.sv = BetaD.sv,
              BetaV.sv = BetaV.sv,
              
              year.p = year.p,
              sd.p = sd.p,
              
              Xi.sv = Xi.sv,
              Epsilon.sv = Epsilon.sv,
              Gamma.sv = Gamma.sv,
              Tau.sv = Tau.sv,
              
              sv = sv,
              b = b,
              s.PY = s.PY,
              s.YAF = s.YAF,
              s.SA = s.SA,
              s.AD = s.AD,
              
              YAF = YAF,
              SA = SA,
              AD = AD,
              Ntot = Ntot
              ))
  
}

# test <- simulateInits(N.year = 40, N.age = 20)
# test

