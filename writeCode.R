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
  # nYAFa = number of young-at-foot of mothers of each age class in the population
  # nSA = number of subadults (1 year old) in the population
  # nAD = number of adults (2 through 19 years old) in the population
  # nTOT = number of female kangaroos from YAF age onwards in the population
  
  # S = year & age-specific survival probabilities from the Cormack-Jolly-Seber model
  # Bt = year-specific breeding rate, or probability of a female of any age producing a jellybean
  # Ra = year & age-specific probability of successfully carrying a jellybean to its 1st Sept as a YAF
  # sYAF = survival of young-at-foot to the 1st of Sept when they are 1 year old (ageC 1 in the CJS model)
  # sSA = survival of 1 year-olds to the 1st of Sept when they are 2 years old (ageC 2 in the CJS model)
  # sAD = survival of adult females of any age from one 1st of Sept to the next (ageC 2, 3, 4 & 5)
  
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
  # Sigma.O = standard deviation of effect of year on probability of observation (was sd.p)
  
  # Mu.B = mean breeding rate, or probability of a female producing a jellybean
  # Mu.R = mean probability of successfully turning a jellybean into a YAF
  
  # BetaD.R = covariate effect of density (D) on reproductive success (R)
  # BetaV.R = covariate effect of vegetation (V) on reproductive success (R)
  # BetaW.R = covariate effect of winter severity (W) on reproductive success (R)
  
  # EpsilonI.R = random effect of mother's identity (I) on reproductive success (Ri)
  # EpsilonT.R = random effect of year (T) on reproductive success (Ri & Ra)
  # EpsilonT.B = random effect of year (T) on breeding rate (Bt)
  
  # SigmaI.R = standard deviation of effect of mother's identity (I) on reproductive success (Ri)
  # SigmaT.R = standard deviation of effect of year (T) on reproductive success (Ri & Ra)
  # SigmaT.B = standard deviation of effect of year (T) on breeding rate (Bt)
  
  # ab = yearly abundance in number of kangaroos, supplied for most years with 1-2 NAs
  # propF = yearly proportion of observations representing females, supplied up to 2019
  
  
  ## Set up --------------------------------------------------------------------
  
  # load packages
  suppressPackageStartupMessages(library(tidyverse))
  library(lubridate)
  library(nimble)
  
  
  ## Model ---------------------------------------------------------------------
  
  myCode = nimbleCode({
    
    
    ## MISSING VALUES
    ## -------------------------------------------------------------------------
    
    if(envEffectsS || envEffectsR){
      for(t in 1:(nYear-1)){
        dens.hat[t] ~ dnorm(dens[t], sd = densE[t])
        veg.hat[t]  ~ dnorm(veg[t], sd = vegE[t])
        win.hat[t]  ~ dnorm(win[t], sd = 1)
      }
      
      for(m in 1:nNoVeg){
        veg[m] ~ dnorm(0, sd = 2)
      }
      
      for(m in 1:nNoDens){
        dens[m] ~ dnorm(0, sd = 2)
      }
      
      for(m in 1:nNoWin){
        win[m] ~ dnorm(0, sd = 2)
      }
      
      for(p in 1:nNoProp){
        propF[p] ~ T(dnorm(0.8, 0.2), 0, 1)
      }
    }
    
    
    ## POPULATION MODEL
    ## -------------------------------------------------------------------------
    
    nAD[1, 1:nYear] <- 0
    nTOT[1] <- nYAF[1] + nSA[1] + sum(nAD[2:nAge, 1])
    
    for(t in 1:(nYear-1)){
      # survival & birthdays
      nSA[t+1] ~ dbin(sYAF[t], nYAF[t])
      nAD[2, t+1] ~ dbin(sSA[t], nSA[t])
      for(a in 3:nAge){
        nAD[a, t+1] ~ dbin(sAD[a-1, t], nAD[a-1, t])
      }
      # then reproduction
      for(a in 3:nAge){
        nYAFa[a, t+1] ~ dbin(0.5 * Bt[t] * Ra[a-1, t], nAD[a-1, t])
      }
      nYAF[t+1] <- sum(nYAFa[3:nAge, t+1]) # number of female YAFs
      nTOT[t+1] <- nYAF[t+1] + nSA[t+1] + sum(nAD[2:nAge, t+1])
    }
    
    # priors
    for(t in 1:(nYear-1)){
      sYAF[t]   <- S[1, t]
      sSA[t]    <- S[2, t]
      sAD[1, t] <- 0 # don't exist
      sAD[2, t] <- S[2, t]
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
    
    
    ## ABUNDANCE MODEL
    ## -------------------------------------------------------------------------
    
    #### Likelihood ####
    for(t in 1:nYear){
      ab[t] ~ dnorm((nTOT[t] / propF[t]), sd = abE[t])
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
        if(envEffectsS){
          logit(S[a, t]) <- BetaA.S[a] +
            BetaD.S[a] * dens.hat[t] +
            BetaV.S[a] * veg.hat[t] +
            # BetaDV.S[a] * (dens.hat[t] * veg.hat[t]) +
            # BetaVR.S[a] * (veg.hat[t] / dens.hat[t]) +
            Gamma.S[t, a]
        }else{
          logit(S[a, t]) <- BetaA.S[a] +
            Gamma.S[t, a]
        }
      }
    }

    # observation function
    for(t in 1:nYear){
      Epsilon.O[t] ~ dnorm(0, sd = Sigma.O)
      logit(O[t]) <- logit(Mu.O) + Epsilon.O[t]
    }
    
    #### Priors ####
    # for fixed effects
    for(a in 1:nAgeC){
      BetaA.S[a] ~ dunif(-5, 5)
      if(envEffectsS){
        BetaD.S[a] ~ dunif(-5, 5)
        BetaV.S[a] ~ dunif(-5, 5)
        # BetaDV.S[a] ~ dunif(-5, 5)
      }
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
    # yearly birth rate
    for(x in 1:nR){
      B[x] ~ dbern(Bi[x])
      logit(Bi[x]) <- logit(Mu.B) +
        EpsilonT.B[year.R[x]]
    }

    for(t in 1:(nYear-1)){
      logit(Bt[t]) <- logit(Mu.B) + EpsilonT.B[t]
    }

    # individual RS function
    for(x in 1:nR){
      if(envEffectsR){
        R[x] ~ dbern(Ri[x])
        logit(Ri[x]) <- logit(Mu.R[age.R[x]]) +
          BetaD.R * dens.hat[year.R[x]] +
          BetaV.R * veg.hat[year.R[x]] +
          BetaW.R * win.hat[year.R[x]] +
          EpsilonI.R[id.R[x]] +
          EpsilonT.R[year.R[x]]
      }else{
        R[x] ~ dbern(Ri[x])
        logit(Ri[x]) <- logit(Mu.R[age.R[x]]) +
          EpsilonI.R[id.R[x]] +
          EpsilonT.R[year.R[x]]
      }
    }

    # age-specific RS function
    # use parameters estimated from individual data above
    # to predict age-specific reproductive success (Ra) here!
    Mu.R[1] <- 0
    for(a in 1:nAge){
      for(t in 1:(nYear-1)){
        if(envEffectsR){
          logit(Ra[a, t]) <- logit(Mu.R[a]) +
            BetaD.R * dens.hat[t] +
            BetaV.R * veg.hat[t] +
            BetaW.R * win.hat[t] +
            EpsilonT.R[t]
        }else{
          logit(Ra[a, t]) <- logit(Mu.R[a]) +
            EpsilonT.R[t]
        }
      }
    }

    ##### Priors ####
    # priors for fixed effects
    for(a in 2:nAge){
      Mu.R[a] ~ dunif(0, 1)
    }
    Mu.B ~ dunif(0, 1)

    if(envEffectsR){
      BetaD.R ~ dunif(-5, 5)
      BetaV.R ~ dunif(-5, 5)
      BetaW.R ~ dunif(-5, 5)
    }
    
    # priors for random effects
    for(i in 1:nID.R){
      EpsilonI.R[i] ~ dnorm(0, sd = SigmaI.R)
    }

    for(t in 1:(nYear-1)){
      EpsilonT.R[t] ~ dnorm(0, sd = SigmaT.R)
      EpsilonT.B[t] ~ dnorm(0, sd = SigmaT.B)
    }

    # priors for sigma
    SigmaI.R <- 0
    # SigmaI.R ~ dunif(0, 100)
    SigmaT.R ~ dunif(0, 100)
    SigmaT.B ~ dunif(0, 100)
    
  }) # nimbleCode
  
}

