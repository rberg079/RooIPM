#' Simulate initial values for IPM
#'
#' @param ntimes integer. Number of time steps in the model. ntimes = 17 by default.
#' @param nADs integer. Number of adult age classes, or maximum age. nADs = 22 by default.
#' @param nAgeC integer. Number of age classes in CJS model. nAge_CJS = 5 by default.
#'
#' @returns list containing all initial values needed for IPM.
#' @export
#'
#' @examples

simulateInits <- function(ntimes = 17, nADs = 22, nAgeC = 5,
                          dens, veg, nNoAge){
  
  
  ## Simulate latent states for input data -------------------------------------
  
  # True environmental conditions (CJS)
  dens.hat <- ifelse(is.na(dens), rnorm(length(dens), 0, .1), dens)
  veg.hat <- ifelse(is.na(veg), rnorm(length(veg), 0, .1), veg)
  
  # Unobserved ages (CJS)
  ageM <- sample(3:8, size = nNoAge, replace = T)
  
  # CJS latent states (mu1, mu2)
  # TODO: simulate these!!
  
  
  ## Simulate vital rate hyperparameters ---------------------------------------
  
  # Covariate effects on survival (CJS)
  B.age <- rnorm(nAgeC, 0, 0.25)
  B.veg <- rnorm(nAgeC, 0, 0.5)
  B.dens <- rnorm(nAgeC, 0, 1)
  # B.densVeg <- rnorm(nAgeC, 0, 1)
  
  
  ## Simulate vital rate random effects ----------------------------------------
  
  # Correlated age-year random effects on survival (CJS)
  # variance-covariance matrix
  zero <- rep(0, nAgeC)
  xi <- rnorm(nAgeC, 1, 0.1)
  
  eps.raw <- matrix(rnorm((ntimes-1)*nAgeC, 0, 0.1),
                    ncol = nAgeC, nrow = (ntimes-1))
  
  gamma <- matrix(NA, ncol = nAgeC, nrow = ntimes-1)
  
  for (t in 1:(ntimes-1)){
    for (a in 1:nAgeC){
      gamma[t,a] <- xi[a] * eps.raw[t,a]
    } # a
  } # t
  
  Tau.raw <- diag(nAgeC) + rnorm(nAgeC^2, 0, 0.1)
  Tau.raw <- inverse((Tau.raw + t(Tau.raw))/2) # should this be Sigma.raw?
  Sigma.raw <- inverse(Tau.raw)
  
  # # uniform covariance matrix
  # zero <- rep(0, nAgeC)
  # sd.yr <- runif(nAgeC, 0, 1)
  # cor.yr <- diag(nAgeC) + 0.01
  # TODO: how to initialize gamma here?

  
  ## Simulate yearly vital rates -----------------------------------------------
  
  # Age class-specific survival rates (CJS)
  s <- matrix(NA, nrow = nAgeC, ncol = ntimes-1)
  
  for(a in 1:nAgeC){
    for(t in 1:(ntimes-1)){
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
  b <- runif(ntimes, 0.5, 1) # raw means span 0.58-0.92
  
  # Survival of PYs
  # to 1st Sept 1 when they become YAFs
  s.PY <- runif(ntimes, 0.1, 1) # raw means span 0.24-0.95
  
  # Survival of YAFs
  # to 2nd Sept 1 when they become SA1s
  # 1st age class considered in our published CJS model
  # s.YAF <- runif(ntimes, 0, 1) # raw means span 0.01-0.88
  s.YAF <- s[1, 1:(ntimes-1)]
  
  # Survival of SA1s to SA2 & SA2 to AD3
  # 2nd age class in our published CJS model
  # s.SA <- matrix(runif(2 * (ntimes), 0.5, 1), nrow = 2)
  s.SA <- rbind(s[2, 1:(ntimes-1)], s[2, 1:(ntimes-1)])
  
  # Survival of all ADs
  s.AD <- matrix(0, nrow = nADs, ncol = ntimes-1)  # 0s for AD1 & AD2
  
  for(a in 3:6){
    s.AD[a, 1:(ntimes-1)] <- s[3, 1:(ntimes-1)]    # prime-age
  }
  
  for(a in 7:9){
    s.AD[a, 1:(ntimes-1)] <- s[4, 1:(ntimes-1)]    # pre-senescent
  }
  
  for(a in 10:nADs){
    s.AD[a, 1:(ntimes-1)] <- s[5, 1:(ntimes-1)]    # senescent
  }
  
  
  ## Simulate observation parameters -------------------------------------------
  
  # Recapture probabilities (CJS)
  mean.p <- runif(1, 0.6, 1)
  year.p <- rnorm(ntimes, 0, 0.2)
  sd.p <- rnorm(1, 0.2, 0.1)
  
  
  ## Simulate initial population sizes -----------------------------------------
  
  # Actual numbers in 2008:
  # 5 female YAFs in Sept, 6 SA1s, 5 SA2s, 21 adults
  # Wendy estimated 22.6% of the population was marked
  
  YAF    <- c(5*5, rep(NA, times = ntimes-1))
  SA     <- matrix(NA, nrow = 2, ncol = ntimes)
  SA[,1] <- c(6*5, 5*5)
  
  AD     <- matrix(NA, nrow = nADs, ncol = ntimes)
  AD[,1] <- c(0, 0, rep(2*5, times = 8), rep(1*5, times = 8))
  
  Ntot   <- c(YAF[1] + sum(SA[1:2,1]) + sum(AD[3:nADs,1]),
              rep(NA, times = ntimes-1))
  
  
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

# test <- simulateInits(ntimes = 40, nADs = 20)
# test

