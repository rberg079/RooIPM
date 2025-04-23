#' Simulate initial values for IPM
#'
#' @param N.year integer. Number of time steps in the model. N.year = 17 by default.
#' @param N.age integer. Number of adult age classes, or maximum age. N.age = 22 by default.
#' @param N.ageC integer. Number of age classes in CJS model. N.ageC = 5 by default.
#' @param dens vector. Where the density data is stored.
#' @param veg vector. Where the vegetation data is stored.
#' @param N.noAge integer. Number of individuals of unknown age in the analysis.
#'
#' @returns list containing all initial values needed for the IPM.
#' @export
#'
#' @examples

simulateInits <- function(N.year = 17, N.age = 22, N.ageC = 5,
                          dens, veg, N.noAge){
  
  
  ## Simulate latent states for input data -------------------------------------
  
  # True environmental conditions (CJS)
  dens.hat <- ifelse(is.na(dens), rnorm(length(dens), 0, .1), dens)
  veg.hat <- ifelse(is.na(veg), rnorm(length(veg), 0, .1), veg)
  
  # Unobserved ages (CJS)
  ageM <- sample(3:8, size = N.noAge, replace = T)
  
  # CJS latent states (mu1, mu2)
  # TODO: simulate these!!
  
  
  ## Simulate vital rate hyperparameters ---------------------------------------
  
  # Covariate effects on survival (CJS)
  B.age <- rnorm(N.ageC, 0, 0.25)
  B.veg <- rnorm(N.ageC, 0, 0.5)
  B.dens <- rnorm(N.ageC, 0, 1)
  # B.densVeg <- rnorm(N.ageC, 0, 1)
  
  
  ## Simulate vital rate random effects ----------------------------------------
  
  # Correlated age-year random effects on survival (CJS)
  # variance-covariance matrix
  zero <- rep(0, N.ageC)
  xi <- rnorm(N.ageC, 1, 0.1)
  
  eps.raw <- matrix(rnorm((N.year-1)*N.ageC, 0, 0.1),
                    ncol = N.ageC, nrow = (N.year-1))
  
  gamma <- matrix(NA, ncol = N.ageC, nrow = N.year-1)
  
  for (t in 1:(N.year-1)){
    for (a in 1:N.ageC){
      gamma[t,a] <- xi[a] * eps.raw[t,a]
    } # a
  } # t
  
  Tau.raw = diag(N.ageC) + rnorm(N.ageC^2, 0, 0.1)
  Tau.raw = matrix(Tau.raw, nrow = N.ageC)
  Tau.raw = (Tau.raw + t(Tau.raw)) / 2
  
  # # uniform covariance matrix
  # zero <- rep(0, N.ageC)
  # sd.yr <- runif(N.ageC, 0, 1)
  # cor.yr <- diag(N.ageC) + 0.01
  # TODO: how to initialize gamma here?

  
  ## Simulate yearly vital rates -----------------------------------------------
  
  # Age class-specific survival rates (CJS)
  s <- matrix(NA, nrow = N.ageC, ncol = N.year-1)
  
  for(a in 1:N.ageC){
    for(t in 1:(N.year-1)){
      s[a, t] <- plogis(
        B.age[a] +
        # B.dens[a] * dens.hat[t] +
        # B.veg[a] * veg.hat[t] +
        # B.densVeg[a] * (dens.hat[t] * veg.hat[t]) +
        # B.vegRoo[a] * (veg.hat[t] / dens.hat[t]) +
        gamma[t, a])
    }
  }

  # Breeding rate
  b <- runif(N.year-1, 0.5, 1) # raw means span 0.58-0.92
  
  # Survival of PYs
  # to 1st Sept 1 when they become YAFs
  s.PY <- runif(N.year-1, 0.1, 1) # raw means span 0.24-0.95
  
  # Survival of YAFs
  # to 2nd Sept 1 when they become SA1s
  # 1st age class considered in our published CJS model
  # s.YAF <- runif(N.year, 0, 1) # raw means span 0.01-0.88
  s.YAF <- s[1, 1:(N.year-1)]
  
  # Survival of SA1s to SA2 & SA2 to AD3
  # 2nd age class in our published CJS model
  # s.SA <- matrix(runif(2 * (N.year), 0.5, 1), nrow = 2)
  s.SA <- rbind(s[2, 1:(N.year-1)], s[2, 1:(N.year-1)])
  
  # Survival of all ADs
  s.AD <- matrix(0, nrow = N.age, ncol = N.year-1)
  
  for(a in 3:6){
    s.AD[a, 1:(N.year-1)] <- s[3, 1:(N.year-1)]
  }
  
  for(a in 7:9){
    s.AD[a, 1:(N.year-1)] <- s[4, 1:(N.year-1)]
  }
  
  for(a in 10:N.age){
    s.AD[a, 1:(N.year-1)] <- s[5, 1:(N.year-1)]
  }
  
  
  ## Simulate observation parameters -------------------------------------------
  
  # Recapture probabilities (CJS)
  mean.p <- runif(1, 0.6, 1)
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
  
  return(list(dens.hat = dens.hat,
              veg.hat = veg.hat,
              ageM = ageM,
              
              B.age = B.age,
              B.dens = B.dens,
              B.veg = B.veg,
              
              mean.p = mean.p,
              year.p = year.p,
              sd.p = sd.p,
              
              xi = xi,
              eps.raw = eps.raw,
              gamma = gamma,
              Tau.raw = Tau.raw,
              Sigma.raw = Sigma.raw,
              
              # cor.yr = cor.yr,
              # sd.yr = sd.yr,
              
              s = s,
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

