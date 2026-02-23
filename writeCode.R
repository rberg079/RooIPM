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
  # nYear = number of years in the population model (was ntimes or N.year)
  # nAge = number of ages or maximum age in the population model
  # nAgeC.S = number of age classes in the survival model (was nAgeC)
  # nAgeC.R = number of age classes in the reproductive success model
  
  # nNoAge = number of individuals for which age is unknown
  # nNoDens = number of years for which population density is unknown
  # nNoVeg = number of years for which available vegetation is unknown
  # nNoWin = number of years for which winter severity is unknown
  
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
  
  # Mu.Sp = mean latent state for CJS model (was mu1)
  # Mu.Op = mean latent observation for CJS model (was mu2)
  
  # Mu.S = age-specific mean probability of survival (S)
  # BetaD.S = covariate effect of density (D) on survival (S) (was B.dens)
  # BetaV.S = covariate effect of vegetation (V) on survival (S) (was B.veg)
  # BetaW.S = covariate effect of weather harshness (W) on survival (S)
  
  # dens.true = "true" yearly population density, from which the observed value was hypothetically sampled
  # dens.cov = centered "true" yearly population density, for its use as a covariate in survival & reproductive success models
  # veg.true = "true" yearly available vegetation, from which the observed value was hypothetically sampled
  # noAge = indexes of individuals who are of unknown age in the survival model
  # ageM = estimated ages of unknown-aged individuals in the survival model
  
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
  # Mu.R = mean probability of successfully turning a jellybean into a YAF, in the reproductive success model
  
  # BetaD.R = covariate effect of density (D) on reproductive success (R)
  # BetaV.R = covariate effect of vegetation (V) on reproductive success (R)
  # BetaW.R = covariate effect of weather harshness (W) on reproductive success (R)
  
  # EpsilonI.R = random effect of mother's identity (I) on reproductive success (Ri)
  # EpsilonT.R = random effect of year (T) on reproductive success (Ri & Ra)
  # EpsilonT.B = random effect of year (T) on breeding rate (Bt)
  
  # XiI.R = latent standard normal scale of the random effect of mother's identity (I) on reproductive success (R)
  # XiT.R = latent standard normal scale of the random effect of year (T) on reproductive success (R)
  # XiT.B = latent standard normal scale of the random effect of year (T) on breeding rate (Bt)
  
  # SigmaI.R = standard deviation of effect of mother's identity (I) on reproductive success (Ri)
  # SigmaT.R = standard deviation of effect of year (T) on reproductive success (Ri & Ra)
  # SigmaT.B = standard deviation of effect of year (T) on breeding rate (Bt)
  
  # propF = yearly proportion of observations representing females, supplied up to 2019
  
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
    
    if(envEffectsS || envEffectsR){
      for(t in 1:(nYear-1)){
        # CRN: This is the data likelihood for the population density estimates. We need it.
        dens[t] ~ dnorm(dens.true[t], sd = densE[t])
        
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
    # initial population sizes
    nYF[1] <- initN.YF
    nSA[1] <- initN.SA
    
    # Poisson priors centered on your simulated/expected starting values
    initN.YF ~ dpois(100) 
    initN.SA ~ dpois(120) 
    
    for(a in 2:nAge){
      nAD[a, 1] <- initN.AD[a]
      initN.AD[a] ~ dpois(80) 
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
    for(i in 1:nID.S){
      for(t in (first[i] + 1):last[i]){
        # state process
        state[i, t] ~ dbern(Mu.Sp[i, t])
        Mu.Sp[i, t] <- S[ageC.S[age.S[i, t-1]], t-1] * state[i, t-1]
        
        # observation process
        obs[i, t] ~ dbern(Mu.Op[i, t])
        Mu.Op[i, t] <- O[t] * state[i, t]
      }
      
      # CRN: I think we are ready to try what we can gain by marginalizing this likelihood. 
      # Try to implement nimbleEcology::dCJS_vv()
      # Documentation here: https://cran.r-project.org/web/packages/nimbleEcology/vignettes/Introduction_to_nimbleEcology.html
    }
    
    #### Constraints ####
    # survival function
    for(a in 1:nAgeC.S){
      for(t in 1:(nYear-1)){
        if(envEffectsS){
          logit(S[a, t]) <- logit(Mu.S[a]) +
            BetaD.S * dens.cov[t] +
            BetaV.S * veg.true[t] +
            BetaW.S * win.true[t] +
            # Gamma.S[t, a] +
            EpsilonT.S[t]
        }else{
          logit(S[a, t]) <- logit(Mu.S[a]) +
            # Gamma.S[t, a] +
            EpsilonT.S[t]
        }
      }
    }

    # observation function
    for(t in 1:nYear){
      EpsilonT.O[t] ~ dnorm(0, sd = SigmaT.O)
      logit(O[t]) <- logit(Mu.O) + EpsilonT.O[t]

      # CRN: The model struggles with the estimation of random effects on O.
      # Can try two alternative models here and see how it affects performance.

      # # 1) Constant model:
      # EpsilonT.O[t] <- 0
    }
    
    # # 2) Autoregressive model
    # EpsilonT.O[1] ~ dnorm(0, sd = SigmaT.O)
    # logit(O[1]) <- logit(Mu.O) + EpsilonT.O[1]
    # 
    # for(t in 2:nYear){
    #   EpsilonT.O[t] ~ dnorm(0, sd = SigmaT.O)
    #   logit(O[t]) <- logit(O[t-1]) + EpsilonT.O[t]
    # }
    
    #### Priors ####
    # survival
    for(a in 1:nAgeC.S){
      Mu.S[a] ~ dunif(0, 1)
    }
    
    # # for age-dependent fixed effects
    # if(envEffectsS){
    #   for(a in 1:nAgeC.S){
    #     BetaD.S[a] ~ dunif(-5, 5)
    #     BetaV.S[a] ~ dunif(-5, 5)
    #     BetaW.S[a] ~ dunif(-5, 5)
    #   }
    # }
    
    # for age-independent fixed effects
    if(envEffectsS){
      BetaD.S ~ dunif(-5, 5)
      BetaV.S ~ dunif(-5, 5)
      BetaW.S ~ dunif(-5, 5)
    }
    
    # # for age-dependent random effects
    # # variance-covariance matrix
    # for(a in 1:nAgeC.S){
    #   zero[a] <- 0
    #   Xi.S[a] ~ dunif(0, 2)
    # }
    # 
    # for(t in 1:(nYear-1)){
    #   Epsilon.S[t, 1:nAgeC.S] ~ dmnorm(zero[1:nAgeC.S], Tau.S[1:nAgeC.S, 1:nAgeC.S])
    #   for(a in 1:nAgeC.S){
    #     Gamma.S[t, a] <- Xi.S[a] * Epsilon.S[t, a]
    #   }
    # }
    # 
    # # precision matrix
    # Tau.S[1:nAgeC.S, 1:nAgeC.S] ~ dwish(W[1:nAgeC.S, 1:nAgeC.S], DF)
    # Sigma.S[1:nAgeC.S, 1:nAgeC.S] <- inverse(Tau.S[1:nAgeC.S, 1:nAgeC.S])
    
    # for age-independent random effects
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
    for(x in 1:nR){
      B[x] ~ dbern(Bt[year.R[x]])
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
          BetaV.R * veg.true[year.R[x]] +
          BetaW.R * win.true[year.R[x]] +
          EpsilonI.R[id.R[x]] +
          EpsilonT.R[year.R[x]]
      }else{
        R[x] ~ dbern(Ri[x])
        logit(Ri[x]) <- logit(Mu.R[ageC.R[age.R[x]]]) +
          EpsilonI.R[id.R[x]] +
          EpsilonT.R[year.R[x]]
      }
    }
    
    # CRN: Posteriors and traceplots indicate that this model is too complex. 
    # We may have to figure out which aspect is pushing it too much. 
    # Try to run some alternatives: 
    # a) No individual random effect
    # b) No density effect
    # c) No covariate effects at all

    # age-specific RS function
    # use parameters estimated from individual data above
    # to predict age-specific reproductive success (Ra) here!
    # Mu.R[1] <- 0
    for(a in 1:nAgeC.R){
      for(t in 1:(nYear-1)){
        if(envEffectsR){
          logit(Ra[a, t]) <- logit(Mu.R[a]) +
            BetaD.R * dens.cov[t] +
            BetaV.R * veg.true[t] +
            BetaW.R * win.true[t] +
            EpsilonT.R[t]
        }else{
          logit(Ra[a, t]) <- logit(Mu.R[a]) +
            EpsilonT.R[t]
        }
      }
    }

    ##### Priors ####
    # priors for fixed effects
    for(a in 1:nAgeC.R){
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
      XiI.R[i] ~ dnorm(0, sd = 1) # latent standard normal
      # XiI.R[i] <- 0
    }

    for(t in 1:(nYear-1)){
      XiT.R[t] ~ dnorm(0, sd = 1) # latent standard normal
      XiT.B[t] ~ dnorm(0, sd = 1) # latent standard normal
    }
    
    EpsilonI.R[1:nID.R]     <- SigmaI.R * XiI.R[1:nID.R]     # actual random effect
    EpsilonT.R[1:(nYear-1)] <- SigmaT.R * XiT.R[1:(nYear-1)] # actual random effect
    EpsilonT.B[1:(nYear-1)] <- SigmaT.B * XiT.B[1:(nYear-1)] # actual random effect
    
    # this way sampler can move Xi & Sigma independently
    # apparently helps avoid strong correlations between variance parameters & effects, improving mixing
    # & apparently analogous to my already non-centered random effects in the survival model block (ref: chatGPT...)
    
    # SigmaI.R <- 0
    SigmaI.R ~ dunif(0, 10) # scale of the random effect
    SigmaT.R ~ dunif(0, 10) # scale of the random effect
    SigmaT.B ~ dunif(0, 10) # scale of the random effect
    
    # NOTES:
    # Survival: interpret Gamma & summarize variation using Sigma (correlated SDs per age class)
    # Reproduction or Survival: interpret Epsilons & summarize variation using Sigmas (uncorrelated SDs)
    
  }) # nimbleCode
  
}

