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
  
  # nNoDens = number of years for which population density is unknown
  # nNoVeg = number of years for which available vegetation is unknown
  # nNoWin = number of years for which weather harshness is unknown
  # nNoProp = number of years for which proportion of females in the population is unknown
  
  # nYF = number of young-at-foot (near pouch exit at 0 years old) in the population (was nYAF)
  # nYFa = number of young-at-foot of mothers of each age in the population (was nYAFa)
  # nSA = number of subadults (1 year old) in the population
  # nAD = number of adults (2 through 19 years old) in the population
  # nTOT = number of female kangaroos from YAF age onwards in the population
  
  # S = year & age-specific survival probabilities from the Cormack-Jolly-Seber (CJS) model
  # Bt = year-specific probability of a female of any age producing a jellybean, from the reproductive success (RS) model
  # Ra = year & age-specific probability of successfully bringing a jellybean to its 1st Sept as a YAF from the RS model
  # sPY = survival of pouch young/jellybean to their 1st Sept as YAF in the population model (taken from the RS model's Ra)
  # sYF = survival of young-at-foot to their 1st Sept as 1-year-old subadults in the population model (taken from the CJS model's S)
  # sSA = survival of 1-year-old subadults to their 1st Sept as 2-year-olds in the population model (taken from the CJS model's S)
  # sAD = survival of adult females from one Sept to the next in the population model (taken from the CJS model's S)
  
  # Mu.S = age-specific mean probability of survival (S)
  # BetaD.S = covariate effect of density (D) on survival (S) (was B.dens)
  # BetaV.S = covariate effect of vegetation (V) on survival (S) (was B.veg)
  # BetaW.S = covariate effect of weather harshness (W) on survival (S)
  
  # dens.true = "true" yearly population density, from which the observed values were hypothetically sampled
  # dens.cov = centered "true" yearly population density, for its use as a covariate in survival & reproductive success models
  # veg.true = "true" yearly available vegetation, from which the observed values were hypothetically sampled
  # win.true = "true" yearly weather harshness, from which the observed values were hypothetically sampled
  
  # Gamma.S = correlated random effect of year on probability of survival (was gamma)
  # Xi.S = scaling factor for how big the random effect variation is per age class (was xi)
  # Epsilon.S = raw random effect sampled from a multivariate normal with precision Tau.S (was eps.raw)
  # Tau.S = precision matrix (inverse of covariance) describing variance & correlation among age classes (was Tau.raw)
  # Sigma.S = covariance matrix (inverse of precision) describing variance & covariance among age classes (was Sigma.raw)
  
  # EpsilonT.S = random effect of year (T) on survival (S), for uncorrelated random effect
  # XiT.S = latent standard normal scale of the random effect, for uncorrelated random effect
  # SigmaT.S = standard deviation of effect of year (T) on survival, for uncorrelated random effect
  
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
  
  # SigmaI.R = standard deviation of effect of mother's identity (I) on reproductive success (Ri)
  # SigmaT.R = standard deviation of effect of year (T) on reproductive success (Ri & Ra)
  # SigmaT.B = standard deviation of effect of year (T) on breeding rate (Bt)
  
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
    
    for(t in 1:(nYear-1)){
      # CRN: This is the data likelihood for the population density estimates. We need it.
      dens[t] ~ dnorm(dens.true[t], sd = densE[t])
    }
    
    if(envEffectsS || envEffectsR){
      for(t in 1:(nYear-1)){
        # CRN: This is an additional level of stochasticity you impose, assuming that vegetation is observed with a known error.
        # This makes sense to include, but perhaps only once everything else works.
        # veg[t]  ~ dnorm(veg.true[t], sd = vegE[t])
        veg[t] <- veg.true[t]
        
        win[t] <- win.true[t]
      }
      
      for(m in 1:nNoVeg){
        veg.true[noVeg[m]] ~ dnorm(0, sd = 1)
      }
      
      for(m in 1:nNoWin){
        win.true[noWin[m]] ~ dnorm(0, sd = 1)
      }
    }
    
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
      for(a in 3:nAge) nYFa[a, t+1] ~ dbin(0.5 * Bt[t] * sPY[a-1, t], nAD[a-1, t])
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
    log.initN.SA ~ dnorm(log(120), sd = 0.5)
    initN.SA <- round(exp(log.initN.SA))
    nSA[1] <- initN.SA
    
    # adults
    for(a in 2:nAge){
      log.initN.AD[a] ~ dnorm(log(80), sd = 0.5)
      initN.AD[a] <- round(exp(log.initN.AD[a]))
      nAD[a,1] <- initN.AD[a]
    }
    
    
    ## POPULATION DENSITY MODEL
    ## -------------------------------------------------------------------------
    
    #### Likelihood ####
    for(t in 1:nYear){
      dens.true[t] <- (nTOT[t] * propF[t]) / area[t]
      dens.cov[t] <- dens.true[t] - densM # center dens for its use as a covariate
    }
    
    
    ## SURVIVAL MODEL (CJS)
    ## -------------------------------------------------------------------------
    
    #### Likelihood ####
    
    # CRN: I think we are ready to try what we can gain by marginalizing this likelihood. 
    # Try to implement nimbleEcology::dCJS_vv()
    # Documentation here: https://cran.r-project.org/web/packages/nimbleEcology/vignettes/Introduction_to_nimbleEcology.html
    
    if(use_dCJS){
      
      ## Marginalized formulation with nimbleEcology::dCJS
      for(i in 1:nID.S){
        for(t in first[i]:(last[i]-1)){
          S_ind[i, t] <- S[ageC.S[age.S[i, t]], t] 
        }
      }

      # Individuals first captured before last-1 (S_ind = vector)
      for(i in 1:(nID.S.switch-1)){
        obs[i, first[i]:last[i]] ~ dCJS_vv(probSurvive = S_ind[i, first[i]:(last[i]-1)], 
                                           probCapture = O[first[i]:last[i]], 
                                           len = last[i] - first[i] + 1)
      }
      
      # Individuals first captured at last-1 (S_ind = scalar)
      for(i in nID.S.switch:nID.S){
        obs[i, first[i]:last[i]] ~ dCJS_sv(probSurvive = S_ind[i, first[i]], 
                                           probCapture = O[first[i]:last[i]], 
                                           len = last[i] - first[i] + 1)
      }
    }else{
      
      ## Latent state formulation
      for(i in 1:nID.S){
        
        # initial state
        state[i, 1] <- 1
        
        for(t in (first[i] + 1):last[i]){
          # state process
          state[i, t] ~ dbern(S[ageC.S[age.S[i, t-1]], t-1] * state[i, t-1])
          
          # observation process
          obs[i, t] ~ dbern(O[t] * state[i, t])
        }
      }
    }
    
    #### Constraints ####
    # survival function
    for(a in 1:nAgeC.S){
      for(t in 1:(nYear-1)){
        if(envEffectsS){
          logit(S[a, t]) <- logit(Mu.S[a]) +
            BetaD.S * dens.cov[t] * dummy[a] +
            BetaV.S * veg.true[t] * dummy[a] +
            BetaW.S * win.true[t] * dummy[a] +
            EpsilonT.S[t]
        }else{
          logit(S[a, t]) <- logit(Mu.S[a]) +
            EpsilonT.S[t]
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
    
    # fixed effects
    if(envEffectsS){
      BetaD.S ~ dunif(-5, 5)
      BetaV.S ~ dunif(-5, 5)
      BetaW.S ~ dunif(-5, 5)
    }
    
    # random effects
    for(t in 1:(nYear-1)){
      XiT.S[t] ~ dnorm(0, sd = 1) # latent standard normal
      EpsilonT.S[t] <- SigmaT.S * XiT.S[t] # actual random effect
    }
    SigmaT.S ~ dunif(0, 10) # scale of the random effect
    
    # observation
    Mu.O ~ dunif(0.01, 0.99) # or dunif(0, 1)
    SigmaT.O ~ dunif(0.01, 10) # or dunif(0, 10)
    
    
    ## REPRODUCTIVE SUCCESS MODEL
    ## -------------------------------------------------------------------------
    
    #### Likelihood & constraints ####
    # yearly birth rate
    for(x in 1:nB){
      B[x] ~ dbern(Bt[year.B[x]])
    }
    
    for(t in 1:(nYear-1)){
      logit(Bt[t]) <- logit(Mu.B) + EpsilonT.B[t]
    }
    
    # individual RS function
    for(x in 1:nR){
      if(envEffectsR){
        R[x] ~ dbern(Ri[x])
        logit(Ri[x]) <- logit(Mu.R[ageC.R[age.R[x]]]) +
          BetaD.R * dens.cov[year.R[x]] +
          EpsilonI.R[id.R[x]] +
          EpsilonT.R[year.R[x]]
      }else{
        R[x] ~ dbern(Ri[x])
        logit(Ri[x]) <- logit(Mu.R[ageC.R[age.R[x]]]) +
          EpsilonI.R[id.R[x]] +
          EpsilonT.R[year.R[x]]
      }
    }
    
    # age-specific RS function
    # use parameters estimated from individual data above
    # to predict age-specific reproductive success (Ra) here!
    for(a in 1:nAgeC.R){
      for(t in 1:(nYear-1)){
        if(envEffectsR){
          logit(Ra[a, t]) <- logit(Mu.R[a]) +
            BetaD.R * dens.cov[t] +
            EpsilonT.R[t]
        }else{
          logit(Ra[a, t]) <- logit(Mu.R[a]) +
            EpsilonT.R[t]
        }
      }
    }
    
    ##### Priors ####
    # fixed effects
    for(a in 1:nAgeC.R){
      Mu.R[a] ~ dunif(0, 1)
    }
    Mu.B ~ dunif(0, 1)
    
    if(envEffectsR){
      BetaD.R ~ dunif(-5, 5)
    }
    
    # random effects
    for(i in 1:nID.R){
      XiI.R[i] ~ dnorm(0, sd = 1) # latent standard normal
    }
    
    for(t in 1:(nYear-1)){
      XiT.R[t] ~ dnorm(0, sd = 1) # latent standard normal
      XiT.B[t] ~ dnorm(0, sd = 1) # latent standard normal
    }
    
    EpsilonI.R[1:nID.R]     <- SigmaI.R * XiI.R[1:nID.R]     # actual random effect
    EpsilonT.R[1:(nYear-1)] <- SigmaT.R * XiT.R[1:(nYear-1)] # actual random effect
    EpsilonT.B[1:(nYear-1)] <- SigmaT.B * XiT.B[1:(nYear-1)] # actual random effect
    
    # NOTES:
    # this way sampler can move Xi & Sigma independently
    # apparently helps avoid strong correlations between variance parameters & effects, improving mixing
    # & apparently analogous to my already non-centered random effects in the survival model block (ref: chatGPT...)
    
    SigmaI.R ~ dunif(0, 10) # scale of the random effect
    SigmaT.R ~ dunif(0, 10) # scale of the random effect
    SigmaT.B ~ dunif(0, 10) # scale of the random effect
    
    # NOTES:
    # Survival: interpret Gamma & summarize variation using Sigma (correlated SDs per age class)
    # Reproduction or Survival: interpret Epsilons & summarize variation using Sigmas (uncorrelated SDs)
    
  }) # nimbleCode
  
}

