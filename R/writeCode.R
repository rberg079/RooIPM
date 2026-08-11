#' Write model code for IPM
#'
#' @returns an R call object specifying the model structure for the kangaroo IPM.
#' @export
#'
#' @examples

writeCode <- function(){
  
  
  ## Parameters ------------------------------------------------------------------
  
  # nYear = number of years in the population model (was ntimes or N.year)
  # nAge = number of ages or maximum age in the population model
  
  # nB = number of events in the birth rate analysis
  # nR = number of events in the survival-of-pouch-young analysis
  # nID.S = number of unique kangaroos in the survival dataset (was nind)
  # nID.R = number of unique kangaroos in the reproductive success dataset (was N.id)

  # nAgeC.S = number of age classes in the survival model (was nAgeC)
  # nAgeC.R = number of age classes in the reproductive success model
  
  # nYF = number of young-at-foot (near pouch exit at 0 years old) in the population (was nYAF)
  # nYFa = number of young-at-foot of mothers of each age in the population (was nYAFa)
  # nSA = number of subadults (1 year old) in the population
  # nAD = number of adults (2 through 19 years old) in the population
  # nTOT = number of female kangaroos from YAF age onwards in the population
  
  # S = year & age-specific survival probabilities from the Cormack-Jolly-Seber (CJS) model
  # Ba = year & age-specific probability of a female producing a jellybean, from the reproductive success (RS) model
  # Ra = year & age-specific probability of successfully bringing a jellybean to its 1st Sept as a YAF from the RS model
  # sPY = survival of pouch young/jellybean to their 1st Sept as YAF in the population model (taken from the RS model's Ra)
  # sYF = survival of young-at-foot to their 1st Sept as 1-year-old subadults in the population model (taken from the CJS model's S)
  # sSA = survival of 1-year-old subadults to their 1st Sept as 2-year-olds in the population model (taken from the CJS model's S)
  # sAD = survival of adult females from one Sept to the next in the population model (taken from the CJS model's S)
  # BR = year & age-specific birth rate in the population model (taken from the RS model's)
  
  # Mu.S = age-specific mean probability of survival (S)
  
  # BetaD.S = covariate effect of density (D) on survival (S) (was B.dens)
  # BetaD.Sy = covariate effect of density on survival of young-at-foot (y, 0 years)
  # BetaD.Sp = covariate effect of density on survival of prime-aged roos (p, 1-9 years)
  # BetaD.Sp = covariate effect of density on survival of old roos (o, 10+ years)
  
  # BetaV.S = covariate effect of vegetation (V) on survival (S) (was B.veg)
  # BetaV.Sy = covariate effect of vegetation on survival of young-at-foot (y, 0 years)
  # BetaV.Sp = covariate effect of vegetation on survival of prime-aged roos (p, 1-9 years)
  # BetaV.Sp = covariate effect of vegetation on survival of old roos (o, 10+ years)
  
  # dens.true = "true" yearly population density, from which the observed values were hypothetically sampled
  # dens.cov = centered "true" yearly population density, for its use as a covariate in survival & reproductive success models
  # veg.true = "true" yearly available vegetation, from which the observed values were hypothetically sampled
  
  # EpsilonT.S = random effect of year (T) on survival (S)
  # EpsilonT.Sy = random effect of year on survival of young-at-foot (y, 0 years)
  # EpsilonT.Sp = random effect of year on survival of prime-aged roos (p, 1-9 years)
  # EpsilonT.So = random effect of year on survival of old roos (o, 10+ years)
  
  # XiT.S = latent standard normal scale of the random effect of year (T) on survival (S)
  # XiT.Sy = latent standard normal scale of the random effect of year on survival of young-at-foot (y, 0 years)
  # XiT.Sp = latent standard normal scale of the random effect of year on survival of prime-aged roos (p, 1-9 years)
  # XiT.So = latent standard normal scale of the random effect of year on survival of old roos (o, 10+ years)
  
  # SigmaT.S = standard deviation of the random effect of year (T) on survival (S)
  # SigmaT.Sy = standard deviation of the random effect of year on survival of young-at-foot (y, 0 years)
  # SigmaT.Sp = standard deviation of the random effect of year on survival of prime-aged roos (p, 1-9 years)
  # SigmaT.So = standard deviation of the random effect of year on survival of old roos (o, 10+ years)
  
  # O = probability of observation each year in the survival model (was p)
  # Mu.O = mean probability of observation in the survival model (was mu.p)
  # EpsilonT.O = random effect of year on prob. of observation in the survival model (was year.p)
  # SigmaT.O = standard deviation of effect of year on probability of observation in the survival model (was sd.p)
  
  # Mu.B = mean breeding rate, or probability of a female producing a jellybean, in the reproductive success model
  # Mu.R = mean probability of successfully bringing a jellybean to pouch exit, in the reproductive success model
  
  # BetaD.R = covariate effect of density (D) on reproductive success (R)
  
  # EpsilonI.R = random effect of mother's identity (I) on reproductive success (Ri)
  # EpsilonT.R = random effect of year (T) on reproductive success (Ri & Ra)
  # EpsilonT.B = random effect of year (T) on breeding rate (Bt)
  
  # XiI.R = latent standard normal scale of the random effect of mother's identity (I) on reproductive success (R)
  # XiT.R = latent standard normal scale of the random effect of year (T) on reproductive success (R)
  # XiT.B = latent standard normal scale of the random effect of year (T) on breeding rate (Bt)
  
  # SigmaI.R = standard deviation of the random effect of mother's identity (I) on reproductive success (Ri)
  # SigmaT.R = standard deviation of the random effect of year (T) on reproductive success (Ri & Ra)
  # SigmaT.B = standard deviation of the random effect of year (T) on breeding rate (Bt)
  
  # initN.YF = initial population size of young-at-foot
  # initN.SA = initial population size of subadults
  # initN.AD = initial population size of adults
  
  
  ## Set up --------------------------------------------------------------------
  
  # load packages
  suppressPackageStartupMessages(library(tidyverse))
  library(lubridate)
  library(nimble)
  
  
  ## Model ---------------------------------------------------------------------
  
  myCode = nimbleCode({
    
    
    ## MISSING VALUES
    ## -------------------------------------------------------------------------
    
    # density data likelihood
    for(t in 1:nYear){
      dens[t] ~ dnorm(dens.true[t], sd = densE[t])
    }
    
    # data imputation for missing vegetation data
    # assuming observation error with known SD
    if(envEffects.S || envEffects.R){
      for(t in 1:(nYear-1)){
        veg[t] ~ dnorm(veg.true[t], sd = vegE[t])
        veg.true[t] ~ dnorm(0, sd = 1)
      }
    }
    
    # data imputation for missing propF data
    for(m in 1:nNoProp){
      propF[noProp[m]] ~ T(dnorm(0.8, sd = 0.2), 0, 1)
    }
    
    
    ## POPULATION MODEL
    ## -------------------------------------------------------------------------
    
    nAD[1, 1:nYear] <- 0 # 1 y/o adults don't exist
    nTOT[1] <- nYF[1] + nSA[1] + sum(nAD[2:nAge, 1])
    
    for(t in 1:(nYear-1)){
      # survival & birthdays
      nSA[t+1] ~ dbin(sYF[t], nYF[t])
      nAD[2, t+1] ~ dbin(sSA[t], nSA[t])
      for(a in 3:nAge) nAD[a, t+1] ~ dbin(sAD[a-1, t], nAD[a-1, t])
      
      # then reproductive success
      for(a in 3:nAge) nYFa[a, t+1] ~ dbin(0.5 * BR[a-1, t] * sPY[a-1, t], nAD[a-1, t])
      nYF[t+1] <- sum(nYFa[3:nAge, t+1]) # total number of female YAFs every year
      nTOT[t+1] <- nYF[t+1] + nSA[t+1] + sum(nAD[2:nAge, t+1])
    }
    
    # survival by age
    # from estimates by age class
    for(t in 1:(nYear-1)){
      sYF[t] <- S[1, t]
      sSA[t] <- S[2, t]
      sAD[1, t] <- 0
      
      if(ageClasses == 6){
        sAD[2, t] <- S[3, t]
        for(a in 3:6) sAD[a, t] <- S[4, t] # prime-aged
        for(a in 7:9) sAD[a, t] <- S[5, t] # pre-senescent
        for(a in 10:nAge) sAD[a, t] <- S[6, t] # senescent
        
      }else if(ageClasses == 12){
        for(a in 2:11) sAD[a, t] <- S[a+1, t] # younger adults
        for(a in 12:nAge) sAD[a, t] <- S[13, t] # greybeards
        
      }else if(ageClasses == 20){
        for(a in 2:19) sAD[a, t] <- S[a+1, t] # adults
      }
    }
    
    # birth rate by age
    # from estimates by age class
    for(t in 1:(nYear-1)){
      BR[1, t] <- 0 # 1 y/os don't reproduce

      if(ageClasses == 6){
        for(a in 2:4) BR[a, t] <- Ba[a-1, t]
        for(a in 5:6) BR[a, t] <- Ba[4, t]
        for(a in 7:10) BR[a, t] <- Ba[5, t]
        for(a in 11:nAge) BR[a, t] <- Ba[6, t]

      }else if(ageClasses == 12){
        for(a in 2:11) BR[a, t] <- Ba[a-1, t]
        for(a in 12:nAge) BR[a, t] <- Ba[11, t]

      }else if(ageClasses == 20){
        for(a in 2:19) BR[a, t] <- Ba[a-1, t]
      }
    }
    
    # reproductive success by age
    # from estimates by age class
    for(t in 1:(nYear-1)){
      sPY[1, t] <- 0 # 1 y/os don't reproduce
      
      if(ageClasses == 6){
        for(a in 2:4) sPY[a, t] <- Ra[a-1, t]
        for(a in 5:6) sPY[a, t] <- Ra[4, t]
        for(a in 7:10) sPY[a, t] <- Ra[5, t]
        for(a in 11:nAge) sPY[a, t] <- Ra[6, t]
        
      }else if(ageClasses == 12){
        for(a in 2:11) sPY[a, t] <- Ra[a-1, t]
        for(a in 12:nAge) sPY[a, t] <- Ra[11, t]
        
      }else if(ageClasses == 20){
        for(a in 2:19) sPY[a, t] <- Ra[a-1, t]
      }
    }
    
    #### Priors ####
    # young-at-foot
    log.initN.YF ~ dnorm(log(100), sd = 0.5)
    initN.YF <- round(exp(log.initN.YF))
    nYF[1] <- initN.YF
    
    # subadults
    log.initN.SA ~ dnorm(log(100), sd = 0.5)
    initN.SA <- round(exp(log.initN.SA))
    nSA[1] <- initN.SA
    
    # adults
    for(a in 2:nAge){
      log.initN.AD[a] ~ dnorm(log(10), sd = 0.5)
      initN.AD[a] <- round(exp(log.initN.AD[a]))
      nAD[a,1] <- initN.AD[a]
    }
    
    
    ## POPULATION DENSITY MODEL
    ## -------------------------------------------------------------------------
    
    #### Likelihood ####
    for(t in 1:nYear){
      dens.true[t] <- (nTOT[t] * propF[t]) / area[t]
      # center dens for its use as a covariate
      dens.cov[t] <- dens.true[t] - densM 
    }
    
    
    ## SURVIVAL MODEL (CJS)
    ## -------------------------------------------------------------------------
    
    #### Likelihood ####
    if(use_dCJS){
      
      ## marginalized formulation with nimbleEcology::dCJS
      for(i in 1:nID.S){
        for(t in first[i]:(last[i]-1)){
          S_ind[i, t] <- S[ageC.S[age.S[i, t]], t] 
        }
      }
      # roos first captured before last-1 (S_ind = vector)
      for(i in 1:(idx.sv-1)){
        obs[i, first[i]:last[i]] ~ dCJS_vv(probSurvive = S_ind[i, first[i]:(last[i]-1)], 
                                           probCapture = O[first[i]:last[i]], 
                                           len = last[i] - first[i] + 1)
      }
      # roos first captured at last-1 (S_ind = scalar)
      for(i in idx.sv:nID.S){
        obs[i, first[i]:last[i]] ~ dCJS_sv(probSurvive = S_ind[i, first[i]], 
                                           probCapture = O[first[i]:last[i]], 
                                           len = last[i] - first[i] + 1)
      }
    }else{
      
      ## latent state formulation
      for(i in 1:nID.S){
        state[i, 1] <- 1 # initial state
        for(t in (first[i] + 1):last[i]){
          state[i, t] ~ dbern(S[ageC.S[age.S[i, t-1]], t-1] * state[i, t-1])
          obs[i, t] ~ dbern(O[t] * state[i, t])
        }
      }
    }
    
    #### Constraints ####
    # survival function
    for(a in 1:nAgeC.S){
      for(t in 1:(nYear-1)){
        
        # with environmental covariates
        
        if(envEffects.S){
          
          if(splitCovs.S == 3){
            
            if(splitREs.S == 3){
              logit(S[a, t]) <- logit(Mu.S[a]) +
                BetaD.Sy * dens.cov[t] * dummy.Sy[a] +
                BetaD.Sp * dens.cov[t] * dummy.Sp[a] +
                BetaD.So * dens.cov[t] * dummy.So[a] +
                BetaV.Sy * veg.true[t] * dummy.Sy[a] +
                BetaV.Sp * veg.true[t] * dummy.Sp[a] +
                BetaV.So * veg.true[t] * dummy.So[a] +
                EpsilonT.Sy[t] * dummy.Sy[a] +
                EpsilonT.Sp[t] * dummy.Sp[a] +
                EpsilonT.So[t] * dummy.So[a]
              
            }else if(splitREs.S == 1){
              logit(S[a, t]) <- logit(Mu.S[a]) +
                BetaD.Sy * dens.cov[t] * dummy.Sy[a] +
                BetaD.Sp * dens.cov[t] * dummy.Sp[a] +
                BetaD.So * dens.cov[t] * dummy.So[a] +
                BetaV.Sy * veg.true[t] * dummy.Sy[a] +
                BetaV.Sp * veg.true[t] * dummy.Sp[a] +
                BetaV.So * veg.true[t] * dummy.So[a] +
                EpsilonT.S[t]
            }
            
          }else if(splitCovs.S == 2){
            
            if(splitREs.S == 3){
              logit(S[a, t]) <- logit(Mu.S[a]) +
                BetaD.Sy * dens.cov[t] * dummy.Sy[a] +
                BetaD.So * dens.cov[t] * dummy.So[a] +
                BetaV.Sy * veg.true[t] * dummy.Sy[a] +
                BetaV.So * veg.true[t] * dummy.So[a] +
                EpsilonT.Sy[t] * dummy.Sy[a] +
                EpsilonT.Sp[t] * dummy.Sp[a] +
                EpsilonT.So[t] * dummy.So[a]
            
            }else if(splitREs.S == 2){
              logit(S[a, t]) <- logit(Mu.S[a]) +
                BetaD.Sy * dens.cov[t] * dummy.Sy[a] +
                BetaD.So * dens.cov[t] * dummy.So[a] +
                BetaV.Sy * veg.true[t] * dummy.Sy[a] +
                BetaV.So * veg.true[t] * dummy.So[a] +
                EpsilonT.Sy[t] * dummy.Sy[a] +
                EpsilonT.So[t] * dummy.So[a]
              
            }else if(splitREs.S == 1){
              logit(S[a, t]) <- logit(Mu.S[a]) +
                BetaD.Sy * dens.cov[t] * dummy.Sy[a] +
                BetaD.So * dens.cov[t] * dummy.So[a] +
                BetaV.Sy * veg.true[t] * dummy.Sy[a] +
                BetaV.So * veg.true[t] * dummy.So[a] +
                EpsilonT.S[t]
            }
            
          }else if(splitCovs.S == 1){
            
            logit(S[a, t]) <- logit(Mu.S[a]) +
              BetaD.S * dens.cov[t] * dummy.S[a] +
              BetaV.S * veg.true[t] * dummy.S[a] +
              EpsilonT.S[t]
            
          }
          
          # without environmental covariates
          
        }else{
          
          if(splitREs.S == 3){
            logit(S[a, t]) <- logit(Mu.S[a]) +
              EpsilonT.Sy[t] * dummy.Sy[a] +
              EpsilonT.Sp[t] * dummy.Sp[a] +
              EpsilonT.So[t] * dummy.So[a]
            
          }else if(splitREs.S == 2){
            logit(S[a, t]) <- logit(Mu.S[a]) +
              EpsilonT.Sy[t] * dummy.Sy[a] +
              EpsilonT.So[t] * dummy.So[a]
            
          }else if(splitREs.S == 1){
            logit(S[a, t]) <- logit(Mu.S[a]) +
              EpsilonT.S[t]
            
          }
        }
      }
    }
    
    # observation function
    for(t in 1:nYear){
      EpsilonT.O[t] ~ dnorm(0, sd = SigmaT.O)
      logit(O[t]) <- logit(Mu.O) + EpsilonT.O[t]
    }
    
    #### Priors ####
    # survival
    for(a in 1:nAgeC.S){
      Mu.S[a] ~ dunif(0, 1)
    }
    
    if(envEffects.S){
      if(splitCovs.S == 3){
        BetaD.Sy ~ dunif(-5, 5)
        BetaD.Sp ~ dunif(-5, 5)
        BetaD.So ~ dunif(-5, 5)
        BetaV.Sy ~ dunif(-5, 5)
        BetaV.Sp ~ dunif(-5, 5)
        BetaV.So ~ dunif(-5, 5)
        
      }else if(splitCovs.S == 2){
        BetaD.Sy ~ dunif(-5, 5)
        BetaD.So ~ dunif(-5, 5)
        BetaV.Sy ~ dunif(-5, 5)
        BetaV.So ~ dunif(-5, 5)
        
      }else if(splitCovs.S == 1){
        BetaD.S ~ dunif(-5, 5)
        BetaV.S ~ dunif(-5, 5)
      }
    }
    
    if(splitREs.S == 3){
      for(t in 1:(nYear-1)){
        XiT.Sy[t] ~ dnorm(0, sd = 1) 
        XiT.Sp[t] ~ dnorm(0, sd = 1) 
        XiT.So[t] ~ dnorm(0, sd = 1)
        EpsilonT.Sy[t] <- SigmaT.Sy * XiT.Sy[t] 
        EpsilonT.Sp[t] <- SigmaT.Sp * XiT.Sp[t] 
        EpsilonT.So[t] <- SigmaT.So * XiT.So[t] 
      }
      SigmaT.Sy ~ dunif(0, 10) 
      SigmaT.Sp ~ dunif(0, 10) 
      SigmaT.So ~ dunif(0, 10) 
      
    }else if(splitREs.S == 2){
      for(t in 1:(nYear-1)){
        XiT.Sy[t] ~ dnorm(0, sd = 1) 
        XiT.So[t] ~ dnorm(0, sd = 1)
        EpsilonT.Sy[t] <- SigmaT.Sy * XiT.Sy[t] 
        EpsilonT.So[t] <- SigmaT.So * XiT.So[t] 
      }
      SigmaT.Sy ~ dunif(0, 10) 
      SigmaT.So ~ dunif(0, 10) 
      
    }else if(splitREs.S == 1){
      for(t in 1:(nYear-1)){
        XiT.S[t] ~ dnorm(0, sd = 1)
        EpsilonT.S[t] <- SigmaT.S * XiT.S[t] 
      }
      SigmaT.S ~ dunif(0, 10) 
    }
    
    # observation
    Mu.O ~ dunif(0.01, 0.99) # or dunif(0, 1)
    SigmaT.O ~ dunif(0.01, 10) # or dunif(0, 10)
    
    
    ## REPRODUCTIVE SUCCESS MODEL
    ## -------------------------------------------------------------------------
    
    #### Likelihood & constraints ####
    # individual birth rate
    for(x in 1:nB){
      B[x] ~ dbern(Bi[x])
      logit(Bi[x]) <- logit(Mu.B[ageC.R[age.B[x]]]) +
        EpsilonT.B[year.B[x]]
    }

    # age-specific birth rate
    for(a in 1:nAgeC.R){
      for(t in 1:(nYear-1)){
        logit(Ba[a, t]) <- logit(Mu.B[a]) +
          EpsilonT.B[t]
      }
    }
    
    # individual reproductive success
    for(x in 1:nR){
      R[x] ~ dbern(Ri[x])
      
      # with environmental covariates
      
      if(envEffects.R){
        
        if(splitCovs.R == 2){
          
          if(splitREs.R == 2){
            logit(Ri[x]) <- logit(Mu.R[ageC.R[age.R[x]]]) +
              BetaD.Rp * dens.cov[year.R[x]] * dummy.Rp[ageC.R[age.R[x]]] +
              BetaD.Ro * dens.cov[year.R[x]] * dummy.Ro[ageC.R[age.R[x]]] +
              EpsilonI.R[id.R[x]] +
              EpsilonT.Rp[year.R[x]] * dummy.Rp[ageC.R[age.R[x]]] +
              EpsilonT.Ro[year.R[x]] * dummy.Ro[ageC.R[age.R[x]]]
            
          }else if(splitREs.R == 1){
            logit(Ri[x]) <- logit(Mu.R[ageC.R[age.R[x]]]) +
              BetaD.Rp * dens.cov[year.R[x]] * dummy.Rp[ageC.R[age.R[x]]] +
              BetaD.Ro * dens.cov[year.R[x]] * dummy.Ro[ageC.R[age.R[x]]] +
              EpsilonI.R[id.R[x]] +
              EpsilonT.R[year.R[x]]
          }
          
        }else if(splitCovs.R == 1){
          
          if(splitREs.R == 2){
            logit(Ri[x]) <- logit(Mu.R[ageC.R[age.R[x]]]) +
              BetaD.R * dens.cov[year.R[x]] * dummy.R[ageC.R[age.R[x]]] +
              EpsilonI.R[id.R[x]] +
              EpsilonT.Rp[year.R[x]] * dummy.Rp[ageC.R[age.R[x]]] +
              EpsilonT.Ro[year.R[x]] * dummy.Ro[ageC.R[age.R[x]]]
            
          }else if(splitREs.R == 1){
            logit(Ri[x]) <- logit(Mu.R[ageC.R[age.R[x]]]) +
              BetaD.R * dens.cov[year.R[x]] * dummy.R[ageC.R[age.R[x]]] +
              EpsilonI.R[id.R[x]] +
              EpsilonT.R[year.R[x]]
          }
        }
        
        # without environmental covariates
        
      }else{
        
        if(splitREs.R == 2){
          logit(Ri[x]) <- logit(Mu.R[ageC.R[age.R[x]]]) +
            EpsilonI.R[id.R[x]] +
            EpsilonT.Rp[year.R[x]] * dummy.Rp[ageC.R[age.R[x]]] +
            EpsilonT.Ro[year.R[x]] * dummy.Ro[ageC.R[age.R[x]]]
          
        }else if(splitREs.R == 1){
          logit(Ri[x]) <- logit(Mu.R[ageC.R[age.R[x]]]) +
            EpsilonI.R[id.R[x]] +
            EpsilonT.R[year.R[x]]
        }
      }
    }
    
    # age-specific reproductive success
    # uses parameters estimated from individual data above
    # to predict age-specific reproductive success (Ra) here
    for(a in 1:nAgeC.R){
      for(t in 1:(nYear-1)){
        
        # with environmental covariates
        
        if(envEffects.R){
          
          if(splitCovs.R == 2){
            
            if(splitREs.R == 2){
              logit(Ra[a, t]) <- logit(Mu.R[a]) +
                BetaD.Rp * dens.cov[t] * dummy.Rp[a] +
                BetaD.Ro * dens.cov[t] * dummy.Ro[a] +
                EpsilonT.Rp[t] * dummy.Rp[a] +
                EpsilonT.Ro[t] * dummy.Ro[a]
              
            }else if(splitREs.R == 1){
              logit(Ra[a, t]) <- logit(Mu.R[a]) +
                BetaD.Rp * dens.cov[t] * dummy.Rp[a] +
                BetaD.Ro * dens.cov[t] * dummy.Ro[a] +
                EpsilonT.R[t]
            }
            
          }else if(splitCovs.R == 1){
            
            if(splitREs.R == 2){
              logit(Ra[a, t]) <- logit(Mu.R[a]) +
                BetaD.R * dens.cov[t] * dummy.R[a] +
                EpsilonT.Rp[t] * dummy.Rp[a] +
                EpsilonT.Ro[t] * dummy.Ro[a]
              
            }else if(splitREs.R == 1){
              logit(Ra[a, t]) <- logit(Mu.R[a]) +
                BetaD.R * dens.cov[t] * dummy.R[a] +
                EpsilonT.R[t]
            }
          }
          
          # without environmental covariates
          
        }else{
          
          if(splitREs.R == 2){
            logit(Ra[a, t]) <- logit(Mu.R[a]) +
              EpsilonT.Rp[t] * dummy.Rp[a] +
              EpsilonT.Ro[t] * dummy.Ro[a]
            
          }else if(splitREs.R == 1){
            logit(Ra[a, t]) <- logit(Mu.R[a]) +
              EpsilonT.R[t]
          }
        }
      }
    }
    
    ##### Priors ####
    # fixed effects
    for(a in 1:nAgeC.R){
      Mu.R[a] ~ dunif(0, 1)
      Mu.B[a] ~ dunif(0, 1)
    }
    
    if(envEffects.R){
      if(splitCovs.R == 2){
        BetaD.Rp ~ dunif(-5, 5)
        BetaD.Ro ~ dunif(-5, 5)
        
      }else if(splitCovs.R == 1){
        BetaD.R ~ dunif(-5, 5)
      }
    }
    
    # random effects
    for(i in 1:nID.R){
      XiI.R[i] ~ dnorm(0, sd = 1)
    }
    EpsilonI.R[1:nID.R] <- SigmaI.R * XiI.R[1:nID.R]
    SigmaI.R ~ dunif(0, 10)
    
    for(t in 1:(nYear-1)){
      XiT.B[t] ~ dnorm(0, sd = 1)
    }
    EpsilonT.B[1:(nYear-1)] <- SigmaT.B * XiT.B[1:(nYear-1)]
    SigmaT.B ~ dunif(0, 10)
    
    if(splitREs.R == 2){
      for(t in 1:(nYear-1)){
        XiT.Rp[t] ~ dnorm(0, sd = 1)
        XiT.Ro[t] ~ dnorm(0, sd = 1)
        EpsilonT.Rp[t] <- SigmaT.Rp * XiT.Rp[t]
        EpsilonT.Ro[t] <- SigmaT.Ro * XiT.Ro[t]
      }
      SigmaT.Rp ~ dunif(0, 10)
      SigmaT.Ro ~ dunif(0, 10)
      
    }else if(splitREs.R == 1){
      for(t in 1:(nYear-1)){
        XiT.R[t] ~ dnorm(0, sd = 1)
        EpsilonT.R[t] <- SigmaT.R * XiT.R[t]
      }
      SigmaT.R ~ dunif(0, 10)
    }
    
  }) # nimbleCode
  
}

