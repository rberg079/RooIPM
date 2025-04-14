#' Simulate initial values for IPM
#'
#' @param ntimes integer. Number of time steps in the model. ntimes = 17 by default.
#' @param nADs integer. Number of adult age classes, or maximum age. nADs = 22 by default.
#'
#' @returns list containing initial values for b, s.PY, s.YAF, s.SA, s.AD, YAF, SA, AD & Ntot.
#' @export
#'
#' @examples

simulateInits <- function(ntimes = 17, nADs = 22, nAge_CJS = 5){
  
  ## Simulate true environmental conditions
  ##---------------------------------------
  veg.hat <- ifelse(is.na(veg), rnorm(length(veg), 0, .1), veg)
  dens.hat <- ifelse(is.na(dens), rnorm(length(dens), 0, .1), dens)
  
  
  ## Simulate vital rate hyperparameters
  ##------------------------------------
  
  # Covariate effects on survival
  B.age <- rnorm(nAge_CJS, 0, 0.25)
  B.veg <- rnorm(nAge_CJS, 0, 0.5)
  B.dens <- rnorm(nAge_CJS, 0, 1)
  #B.densVeg <- rnorm(nAge_CJS, 0, 1)
           
  
  ## Simulate vital rate random effects
  
  # Correlated age-year REs on survival (CJS)
  
  # Variance-Covariance matrix
  xi <- rnorm(nAge_CJS, 1, 0.1)
  zero <- rep(0, nAge_CJS)
  
  eps.raw <- matrix(rnorm((ntimes-1)*nAge_CJS, 0, 0.1),
                    ncol = nAge_CJS, nrow = (ntimes-1))
  
  # cor.yr <- diag(nAge_CJS)+0.01
  # sd.yr <- runif(nAge_CJS, 0,1)
  
  for (t in 1:(ntimes-1)){
    for (a in 1:nAge_CJS){
      gamma[t,a] <- xi[a] * eps.raw[t,a]
    } #a
  } #t
  
  # Priors for precision matrix
  Tau.raw <- diag(nAge_CJS) + rnorm(nAge_CJS^2,0,0.1)
  Tau.raw <- inverse((Tau.raw + t(Tau.raw))/2)
  Sigma.raw <- inverse(Tau.raw)

  
  
  ## Simulate yearly vital rates ------------------------------------------------------
  ##----------------------------
  
  # Age class-specific survival rates (CJS model)
  s <- matrix(NA, nrow = nAge_CJS, ncol = ntimes-1)
  
  for(a in 1:nAge_CJS){
    for(t in 1:(ntimes-1)){
      logit(s[a,t]) <- B.age[a] +
        dens.hat[t]*B.dens[a] +
        veg.hat[t]*B.veg[a] +
        # (dens.hat[t]*veg.hat[t])*B.densVeg[a] +
        # (veg.hat[t]/dens.hat[t])*B.vegRoo[a] +
        gamma[t,a]
    }
  }

  
  # breeding rate
  # b.exp <- mean(rs$Repro, na.rm = T) # spans 0.58-0.92
  # b <- runif(1, b.exp*0.75, ifelse(b.exp*1.25 < 1, b.exp*1.25, 1))
  b <- runif(ntimes, 0.5, 1)
  
  # survival of PYs
  # to 1st Sept 1 when they become YAFs
  # sPY.exp <- mean(rs$SurvSep1, na.rm = T) # spans 0.24-0.95
  # s.PY <- runif(1, sPY.exp*0.5, ifelse(sPY.exp*2 < 1, sPY.exp*2, 1))
  s.PY <- runif(ntimes, 0.1, 1)
  
  # survival of YAFs
  # to 2nd Sept 1 when they become SA1s
  # 1st age class considered in our published CJS model
  # sYAF.exp <- mean(rs$SurvSep2, na.rm = T) # spans 0.01-0.88
  #s.YAF <- runif(ntimes, 0, 1)
  s.YAF <- s[1, 1:(ntimes-1)]
  
  # survival of SA1s to SA2 & SA2 to AD3
  # 2nd age class in our published CJS model
  #s.SA <- matrix(runif(2 * (ntimes), 0.5, 1), nrow = 2)
  s.SA <- rbind(s[2, 1:(ntimes-1)], s[2, 1:(ntimes-1)])
  
  # survival of all ADs
  s.AD            <- matrix(0, nrow = nADs, ncol = ntimes-1)          # 0s for AD1 & AD2
  
  for(a in 3:6){
    s.AD[a, 1:(ntimes-1)]     <- s[3, 1:(ntimes-1)]  # Prime-age
  }
  
  for(a in 7:9){
    s.AD[7:9, 1:(ntimes-1)]   <- s[4, 1:(ntimes-1)]  # Pre-senescent
  }
  
  for(a in 10:nADs){
    s.AD[10:nADs, ] <- s[5, 1:(ntimes-1)]
  }
  
  
  ## Simulate observation parameters
  ##--------------------------------
  
  # Recapture probabilities
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
  
  
  ## Simulate latent states for input data
  
  # Unobserved ages (CJS)
  ageM <- sample(3:8, size = nNoAge, replace = T)
  
  # CJS latent states (mu1, mu2)
  # --> Not simulated so far. We may have to later. 
  
  ## Assemble myinits list -----------------------------------------------------
  
  return(list(YAF = YAF,
              SA = SA,
              AD = AD,
              Ntot = Ntot,
              b = b,
              s.PY = s.PY,
              s.YAF = s.YAF,
              s.SA = s.SA,
              s.AD = s.AD))
  
}

# test <- simulateInits(ntimes = 40, nADs = 20)
# test

