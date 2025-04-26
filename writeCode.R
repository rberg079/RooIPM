#' Write model code for IPM
#'
#' @returns an R call object specifying the model structure for the kangaroo IPM.
#' @export
#'
#' @examples

writeCode <- function(){
  
  
  ## Parameters ------------------------------------------------------------------
  
  # n = number of events in the reproductive success dataset
  # nIDs = number of unique kangaroos in the survival dataset
  # nIDr = number of unique kangaroos in the reproductive success dataset
  # nYear = number of years in the dataset
  # nAge = number of ages in the analysis (3 through 19 years old, so 17 ages)
  # nAgeC = number of age classes in the analysis (not used so far in RS analysis)
  
  # nNoAge = number of individuals for which age is unknown
  # nNoDens = number of years for which population density is unknown
  # nNoVeg = number of years for which available vegetation is unknown
  # nNoWin = number of years for which winter severity is unknown
  
  # nYAF = number of young-at-foot (near pouch exit at 0 years old) in the population
  # nSA = number of subadults (either 1 or 2 years old) in the population
  # nAD = number of adults (3 through 20 years old) in the population
  # nTOT = number of female kangaroos from YAF age onwards in the population
  
  # b = breeding rate, or proportion of females producing a jellybean
  # svPY = survival of jellybeans to pouch exit (or to the 1st of Sept around that time)
  # svYAF = survival of young-at-foot to the 1st of Sept when they are 1 year old
  # svSA = survival of 1 or 2 year-olds to the following 1st of Sept when they are 2 or 3
  # svAD = survival of adult females of any given age from one 1st of Sept to the next
  # sv = survival of each age class considered in the Cormack-Jolly-Seber model
  
  # Mu.sp = mean latent state for CJS model (was mu1)
  # Mu.op = mean latent observation for CJS model (was mu2)
  
  # BetaA.sv = covariate effect of age (A) on survival (sv) (was B.age)
  # BetaD.sv = covariate effect of density (D) on survival (sv) (was B.dens)
  # BetaV.sv = covariate effect of vegetation (V) on survival (sv) (was B.veg)
  # BetaDV.sv = covariate effect of interacting density & vegetation (DV) on survival (sv) (was B.densVeg)
  # BetaVR.sv = covariate effect of vegetation per capita, or kangaroo (VR) on survival (sv) (was B.vegRoo)
  
  # dens.hat = "true" yearly population density, from which the observed value was hypothetically sampled
  # veg.hat = "true" yearly available vegetation, from which the observed value was hypothetically sampled
  # noAge = indexes of individuals who are of unknown age
  # ageM = estimated ages of unknown-aged individuals
  
  # Gamma.sv = correlated random effect of year on probability of survival (was gamma)
  # Xi.sv = scaling factor for how big the random effect variation is per age class (was xi)
  # Epsilon.sv = raw random effect sampled from a multivariate normal with precision Tau.sv (was eps.raw)
  # Tau.sv = precision matrix (inverse of covariance) describing variance & correlation among age classes (was Tau.raw)
  # Sigma.sv = covariance matrix (inverse of precision) describing variance & covariance among age classes (was Sigma.raw)
  
  # ob = probability of observation in each year considered in the Cormack-Jolly-Seber model (was p)
  # Mu.ob = mean probability of observation (was mu.p)
  # Epsilon.ob = random effect of year on prob. of observation (was year.p)
  # Sigma.ob = standard deviation of effect of year on prob. of observation (was sd.p)
  
  
  ## Set up --------------------------------------------------------------------
  
  # load packages
  library(tidyverse)
  library(lubridate)
  library(nimble)
  
  
  ## Model ---------------------------------------------------------------------
  
  myCode = nimbleCode({
    
    ## POPULATION MODEL
    ## -------------------------------------------------------------------------
    
    nAD[1:2, 1:nYear] <- 0
    nTOT[1] <- nYAF[1] + sum(nSA[1:2, 1]) + sum(nAD[3:(nAge+2), 1])
    
    for(t in 1:(nYear-1)){
      nYAF[t+1] ~ dbin(b[t] * svPY[t], sum(nAD[3:(nAge+2), t])) 
      
      nSA[1, t+1] ~ dbin(svYAF[t], nYAF[t])
      nSA[2, t+1] ~ dbin(svSA[1, t], nSA[1, t])
      
      nAD[3, t+1] ~ dbin(svSA[2, t], nSA[2, t])
      
      for(a in 4:(nAge+2)){
        nAD[a, t+1] ~ dbin(svAD[a, t], nAD[a, t])
      }
      nTOT[t+1] <- nYAF[t+1] + sum(nSA[1:2, t+1]) + sum(nAD[3:(nAge+2), t+1])
    }
    
    # priors
    svAD[1:2, 1:(nYear-1)] <- 0
    
    for(t in 1:(nYear-1)){
      b[t]       ~ dunif(0.5, 1)
      svPY[t]    ~ dunif(0.1, 1)
      svYAF[t]   <- sv[1, t]
      svSA[1, t] <- sv[2, t]
      svSA[2, t] <- sv[2, t]
      
      for(a in 3:6){ # prime-aged
        svAD[a, t] <- sv[3, t]
      }
      
      for(a in 7:9){ # pre-senescent
        svAD[a, t] <- sv[4, t]
      }
      
      for(a in 10:(nAge+2)){ # senescent
        svAD[a, t] <- sv[5, t]
      }
    }
    
    ## CAPTURE-MARK-RECAPTURE MODEL (CJS)
    ## -------------------------------------------------------------------------
    
    #### Likelihood ####
    for(i in 1:nIDs){
      for(t in (first[i] + 1):last[i]){
        # state process
        state[i, t] ~ dbern(Mu.sp[i, t])
        Mu.sp[i, t] <- sv[ageC[age[i, t-1]], t-1] * state[i, t-1]
        
        # observation process
        obs[i, t] ~ dbern(Mu.op[i, t])
        Mu.op[i, t] <- ob[t] * state[i, t]
      }
    }
    
    #### Constraints ####
    # survival function
    for(a in 1:nAgeC){                               
      for(t in 1:(nYear-1)){
        logit(sv[a, t]) <- BetaA.sv[a] +
          BetaD.sv[a] * dens.hat[t] +
          BetaV.sv[a] * veg.hat[t] +
          # BetaDV.sv[a] * (dens.hat[t] * veg.hat[t]) +
          # BetaVR.sv[a] * (veg.hat[t] / dens.hat[t]) +
          Gamma.sv[t, a]
      }
    }
    
    for(t in 1:(nYear-1)){
      dens.hat[t] ~ dnorm(dens[t], sd = densE[t])
      veg.hat[t] ~ dnorm(veg[t], sd = vegE[t])
    }
    
    for(m in 1:nNoVeg){
      # veg[m] <- 0
      vegM[m] ~ dnorm(0, sd = 2)
    }
    
    # estimate missing ages
    for(i in 1:nNoAge){
      ageM[i] ~ T(dnegbin(0.25,1.6), 3, 20)
      age[noAge[i], first[noAge[i]]] <- round(ageM[i]) + 1
      for(t in (first[noAge[i]]+1):nYear){
        age[noAge[i], t] <- age[noAge[i], t-1] + 1
      }
    }
    
    # observation function
    for(t in 1:nYear){
      logit(ob[t]) <- logit(Mu.ob) + Epsilon.ob[t]
      Epsilon.ob[t] ~ dnorm(0, sd = Sigma.ob)
    }
    
    #### Priors ####
    # for fixed effects
    for(a in 1:nAgeC){
      BetaA.sv[a] ~ dnorm(0, sd = 2)
      BetaD.sv[a] ~ dnorm(0, sd = 2)
      BetaV.sv[a] ~ dnorm(0, sd = 2)
      # BetaDV.sv[a] ~ dnorm(0, sd = 2)
    }
    
    # for random effects
    # variance-covariance matrix
    for(i in 1:nAgeC){
      zero[i] <- 0
      Xi.sv[i] ~ dunif(0, 2)
    }
    
    for(t in 1:(nYear-1)){
      Epsilon.sv[t, 1:nAgeC] ~ dmnorm(zero[1:nAgeC], Tau.sv[1:nAgeC, 1:nAgeC])
      for(i in 1:nAgeC){
        Gamma.sv[t, i] <- Xi.sv[i] * Epsilon.sv[t, i]
      }
    }
    
    # precision matrix
    Tau.sv[1:nAgeC, 1:nAgeC] ~ dwish(W[1:nAgeC, 1:nAgeC], DF)
    Sigma.sv[1:nAgeC, 1:nAgeC] <- inverse(Tau.sv[1:nAgeC, 1:nAgeC])
    
    # observation
    Mu.ob ~ dunif(0.01, 0.99) # or dunif(0, 1)
    Sigma.ob ~ dunif(0.01, 10) # or dunif(0, 10)
    
  })
  
}

