#' Write model code for IPM
#'
#' @returns an R call object specifying the model structure for the kangaroo IPM.
#' @export
#'
#' @examples

writeCode <- function(){
  
  
  ## Parameters ------------------------------------------------------------------
  
  # nR = number of events in the reproductive success dataset (was N)
  # nID.S = number of unique kangaroos in the survival dataset (was nind)
  # nID.R = number of unique kangaroos in the reproductive success dataset (was N.id)
  # nYear = number of years in the dataset (was ntimes or N.year)
  # nAge = number of ages or maximum age in the analysis
  # nAgeC = number of age classes in the analysis
  
  # nNoAge = number of individuals for which age is unknown
  # nNoDens = number of years for which population density is unknown
  # nNoVeg = number of years for which available vegetation is unknown
  # nNoWin = number of years for which winter severity is unknown
  
  # nYAF = number of young-at-foot (near pouch exit at 0 years old) in the population
  # nSA = number of subadults (either 1 or 2 years old) in the population
  # nAD = number of adults (3 through 20 years old) in the population
  # nTOT = number of female kangaroos from YAF age onwards in the population
  
  # B = breeding rate, or proportion of females producing a jellybean
  # sPY = survival of jellybeans to pouch exit (or to the 1st of Sept around that time)
  # sYAF = survival of young-at-foot to the 1st of Sept when they are 1 year old
  # sSA = survival of 1 or 2 year-olds to the following 1st of Sept when they are 2 or 3
  # sAD = survival of adult females of any given age from one 1st of Sept to the next
  # S = survival of each age class considered in the Cormack-Jolly-Seber model
  
  # Mu.Sp = mean latent state for CJS model (was mu1)
  # Mu.Op = mean latent observation for CJS model (was mu2)
  
  # BetaA.S = covariate effect of age (A) on survival (s) (was B.age)
  # BetaD.S = covariate effect of density (D) on survival (s) (was B.dens)
  # BetaV.S = covariate effect of vegetation (V) on survival (s) (was B.veg)
  # BetaDV.S = covariate effect of interacting density & vegetation (DV) on survival (s) (was B.densVeg)
  # BetaVR.S = covariate effect of vegetation per capita, or kangaroo (VR) on survival (s) (was B.vegRoo)
  
  # dens.hat = "true" yearly population density, from which the observed value was hypothetically sampled
  # veg.hat = "true" yearly available vegetation, from which the observed value was hypothetically sampled
  # noAge = indexes of individuals who are of unknown age
  # ageM = estimated ages of unknown-aged individuals
  
  # Gamma.S = correlated random effect of year on probability of survival (was gamma)
  # Xi.S = scaling factor for how big the random effect variation is per age class (was xi)
  # Epsilon.S = raw random effect sampled from a multivariate normal with precision Tau.S (was eps.raw)
  # Tau.S = precision matrix (inverse of covariance) describing variance & correlation among age classes (was Tau.raw)
  # Sigma.S = covariance matrix (inverse of precision) describing variance & covariance among age classes (was Sigma.raw)
  
  # O = probability of observation in each year considered in the Cormack-Jolly-Seber model (was p)
  # Mu.O = mean probability of observation (was mu.p)
  # Epsilon.O = random effect of year on prob. of observation (was year.p)
  # Sigma.O = standard deviation of effect of year on prob. of observation (was sd.p)
  
  
  ## Set up --------------------------------------------------------------------
  
  # load packages
  library(tidyverse)
  library(lubridate)
  library(nimble)
  
  # # check that all years are represented in RS data
  # if(setequal(1:17, unique(myData$year.R)) == FALSE){
  #   stop("Some years not represented in rsData.")
  # }
  
  
  ## Model ---------------------------------------------------------------------
  
  myCode = nimbleCode({
    
    ## POPULATION MODEL
    ## -------------------------------------------------------------------------
    
    nAD[1:2, 1:nYear] <- 0
    nTOT[1] <- nYAF[1] + sum(nSA[1:2, 1]) + sum(nAD[3:nAge, 1])
    
    for(t in 1:(nYear-1)){
      # survival & birthdays
      nSA[1, t+1] ~ dbin(sYAF[t], nYAF[t])
      nSA[2, t+1] ~ dbin(sSA[1, t], nSA[1, t])
      nAD[3, t+1] ~ dbin(sSA[2, t], nSA[2, t])
      for(a in 4:nAge){
        nAD[a, t+1] ~ dbin(sAD[a-1, t], nAD[a-1, t])
      }
      # then reproduction
      for(a in 3:nAge){
        nYAFa[a, t+1] ~ dbin(B[t] * Ra[a-1, t], nAD[a-1, t])
      }
      nYAF[t+1] <- sum(nYAFa[3:nAge, t+1])
      nTOT[t+1] <- nYAF[t+1] + sum(nSA[1:2, t+1]) + sum(nAD[3:nAge, t+1])
    }
    
    # priors
    sAD[1:2, 1:(nYear-1)] <- 0
    
    for(t in 1:(nYear-1)){
      B[t]       ~ dunif(0.5, 1)
      sYAF[t]   <- S[1, t]
      sSA[1, t] <- S[2, t]
      sSA[2, t] <- S[2, t]
      
      for(a in 3:6){ # prime-aged
        sAD[a, t] <- S[3, t]
      }
      
      for(a in 7:9){ # pre-senescent
        sAD[a, t] <- S[4, t]
      }
      
      for(a in 10:nAge){ # senescent
        sAD[a, t] <- S[5, t]
      }
    }
    
    ## SURVIVAL MODEL (CJS)
    ## -------------------------------------------------------------------------
    
    #### Likelihood ####
    for(i in 1:nID.S){
      for(t in (first[i] + 1):last[i]){
        # state process
        state[i, t] ~ dbern(Mu.Sp[i, t])
        Mu.Sp[i, t] <- S[ageC[age.S[i, t-1]], t-1] * state[i, t-1]
        
        # observation process
        obs[i, t] ~ dbern(Mu.Op[i, t])
        Mu.Op[i, t] <- O[t] * state[i, t]
      }
    }
    
    #### Constraints ####
    # survival function
    for(a in 1:nAgeC){                               
      for(t in 1:(nYear-1)){
        logit(S[a, t]) <- BetaA.S[a] +
          BetaD.S[a] * dens.hat[t] +
          BetaV.S[a] * veg.hat[t] +
          # BetaDV.S[a] * (dens.hat[t] * veg.hat[t]) +
          # BetaVR.S[a] * (veg.hat[t] / dens.hat[t]) +
          Gamma.S[t, a]
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
      age.S[noAge[i], first[noAge[i]]] <- round(ageM[i]) + 1
      for(t in (first[noAge[i]]+1):nYear){
        age.S[noAge[i], t] <- age.S[noAge[i], t-1] + 1
      }
    }
    
    # observation function
    for(t in 1:nYear){
      logit(O[t]) <- logit(Mu.O) + Epsilon.O[t]
      Epsilon.O[t] ~ dnorm(0, sd = Sigma.O)
    }
    
    #### Priors ####
    # for fixed effects
    for(a in 1:nAgeC){
      BetaA.S[a] ~ dnorm(0, sd = 2)
      BetaD.S[a] ~ dnorm(0, sd = 2)
      BetaV.S[a] ~ dnorm(0, sd = 2)
      # BetaDV.S[a] ~ dnorm(0, sd = 2)
    }
    
    # for random effects
    # variance-covariance matrix
    for(i in 1:nAgeC){
      zero[i] <- 0
      Xi.S[i] ~ dunif(0, 2)
    }
    
    for(t in 1:(nYear-1)){
      Epsilon.S[t, 1:nAgeC] ~ dmnorm(zero[1:nAgeC], Tau.S[1:nAgeC, 1:nAgeC])
      for(i in 1:nAgeC){
        Gamma.S[t, i] <- Xi.S[i] * Epsilon.S[t, i]
      }
    }
    
    # precision matrix
    Tau.S[1:nAgeC, 1:nAgeC] ~ dwish(W[1:nAgeC, 1:nAgeC], DF)
    Sigma.S[1:nAgeC, 1:nAgeC] <- inverse(Tau.S[1:nAgeC, 1:nAgeC])
    
    # observation
    Mu.O ~ dunif(0.01, 0.99) # or dunif(0, 1)
    Sigma.O ~ dunif(0.01, 10) # or dunif(0, 10)
    
    
    ## REPRODUCTIVE SUCCESS MODEL
    ## -------------------------------------------------------------------------
    
    #### Likelihood & constraints ####
    # individual RS function
    Mu.Ri[1:2] <- 0
    Mu.Ra[1:2] <- 0
    
    for(x in 1:nR){
      R[x] ~ dbern(Ri[x])
      logit(Ri[x]) <- logit(Mu.Ri[age.R[x]]) +
        # BetaD.R * dens[year.R[x]] +
        # BetaV.R * veg[year.R[x]] +
        # BetaW.R * win[year.R[x]] +
        EpsilonI.Ri[id.R[x]] +
        EpsilonT.Ri[year.R[x]]
    }
    
    # age-specific RS function
    # use parameters estimated from individual data above
    # to predict age-specific reproductive success (Ra) here!
    for(a in 1:nAge){
      for(t in 1:(nYear-1)){
        logit(Ra[a, t]) <- logit(Mu.Ra[a]) + # Ra used in Pop model
          # BetaD.R * dens[t] +
          # BetaV.R * veg[t] +
          # BetaW.R * win[t] +
          EpsilonT.Ra[t]
      }
    }
    
    ##### Priors ####
    # priors for fixed effects
    for(a in 3:nAge){
      Mu.Ri[a] ~ dunif(0, 1)
      Mu.Ra[a] ~ dunif(0, 1)
    }
    
    # BetaD.R ~ dunif(-2, 2) # could be dunif(-5, 5) if need be
    # BetaV.R  ~ dunif(-2, 2) # could be dunif(-5, 5) if need be
    # BetaW.R  ~ dunif(-2, 2) # could be dunif(-5, 5) if need be
    
    # priors for random effects
    for(i in 1:nID.R){
      EpsilonI.Ri[i] ~ dnorm(0, sd = SigmaI.Ri)
    }
    
    for(t in 1:(nYear-1)){
      EpsilonT.Ri[t] ~ dnorm(0, sd = SigmaT.Ri)
      EpsilonT.Ra[t] ~ dnorm(0, sd = SigmaT.Ra)
    }
    
    # priors for sigma
    SigmaI.Ri ~ dunif(0, 100)
    SigmaT.Ri ~ dunif(0, 100)
    SigmaT.Ra ~ dunif(0, 100)
    
  }) # nimbleCode
  
}

