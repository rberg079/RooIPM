#' Extract posterior samples for vital rates and population-level quantities
#'
#' @param MCMCsamples MCMClist object. The IPM output as returned by nimbleMCMC().
#' @param nYear integer. Number of time steps in the model. nYear = 17 by default.
#' @param nAge integer. Number of ages, or maximum age, in the model. nAge = 17 by default.
#' @param saveList logical. If TRUE, saves the created list object as an RDS to the /results folder.
#' @param testRun logical. If TRUE, will only use the first 100 samples, for computational efficiency.
#'
#' @returns a list of lists containing posterior samples for all vital rates & population-level quantities.
#' The sublist "t" contains time-specific parameters while the sublist "t.mean" contains time-averaged parameters.
#' The latter are needed for evaluating transient sensitivities. 
#' @export
#'
#' @examples
#' 

extractParamSamples <- function(MCMCsamples, nYear = 17, nAge = 18,
                                saveList = FALSE, testRun = FALSE){
  
  # # for testing purposes
  # MCMCsamples <- readRDS('results/IPM_CJSen_RSen_AB_DynDens_dCJS_12_noW_stochV.rds')
  # nYear <- 17
  # nAge <- 18
  # testRun <- FALSE

  ## Set up --------------------------------------------------------------------
  
  # load packages
  suppressPackageStartupMessages(library(tidyverse))
  
  # convert MCMC samples to matrix
  out.mat <- do.call(rbind, lapply(MCMCsamples, as.matrix))
  
  if(testRun){
    out.mat <- out.mat[1:100,]
  }
  
  # extract number of samples
  nSamples <- dim(out.mat)[1]
  
  # prepare arrays
  # time-varying vital rates
  Bt <- matrix(NA, nrow = nSamples, ncol = nYear-1)
  sPY <- array(NA, dim = c(nSamples, nAge, nYear-1))
  sYF <- matrix(NA, nrow = nSamples, ncol = nYear-1)
  sSA <- matrix(NA, nrow = nSamples, ncol = nYear-1)
  sAD <- array(NA, dim = c(nSamples, nAge, nYear-1))
  
  # time-varying population sizes
  nYF <- matrix(NA, nrow = nSamples, ncol = nYear)
  nSA <- matrix(NA, nrow = nSamples, ncol = nYear)
  nAD <- array(NA, dim = c(nSamples, nAge, nYear))
  nTOT <- matrix(NA, nrow = nSamples, ncol = nYear)
  pYF <- matrix(NA, nrow = nSamples, ncol = nYear)
  pSA <- matrix(NA, nrow = nSamples, ncol = nYear)
  pAD <- array(NA, dim = c(nSamples, nAge, nYear))
  lambda <- matrix(NA, nrow = nSamples, ncol = nYear-1)
  
  
  ## Fill samples into arrays --------------------------------------------------
  
  for(i in 1:nSamples){
    for(t in 1:nYear){
      if(t < nYear){
        
        # time-varying vital rates
        Bt[i, t] <- out.mat[i, paste0("Bt[", t, "]")]
        sYF[i, t] <- out.mat[i, paste0("sYF[", t, "]")]
        sSA[i, t] <- out.mat[i, paste0("sSA[", t, "]")]
        
        for(a in 1:nAge){
          sPY[i, a, t] <- out.mat[i, paste0("sPY[", a, ", ", t, "]")]
          sAD[i, a, t] <- out.mat[i, paste0("sAD[", a, ", ", t, "]")]
        }
      }
      
      # time-varying population sizes
      nYF[i, t] <- out.mat[i, paste0("nYF[", t, "]")]
      nSA[i, t] <- out.mat[i, paste0("nSA[", t, "]")]
      
      for(a in 1:nAge){
        nAD[i, a, t] <- out.mat[i, paste0("nAD[", a, ", ", t, "]")]
      }
      
      nTOT[i, t] <- out.mat[i, paste0("nTOT[", t, "]")]
      
      # time-varying proportions of the population
      # represented by each age class
      pYF[i, t] <- out.mat[i, paste0("nYF[", t, "]")] / out.mat[i, paste0("nTOT[", t, "]")]
      pSA[i, t] <- out.mat[i, paste0("nSA[", t, "]")] / out.mat[i, paste0("nTOT[", t, "]")]
      
      for(a in 1:nAge){
        pAD[i, a, t] <- out.mat[i, paste0("nAD[", a, ", ", t, "]")] / out.mat[i, paste0("nTOT[", t, "]")]
      }
      
      # population growth rate
      if(t < nYear){
        lambda[i, t] <- out.mat[i, paste0("nTOT[", t+1, "]")] / out.mat[i, paste0("nTOT[", t, "]")]
      }
    }
  }
  
  
  ## Calculate time-averages ---------------------------------------------------
  
  # vital rates
  Bt.mean <- rowMeans(Bt[, 1:(nYear-1)], na.rm = T)
  sPY.mean <- apply(sPY[, , 1:(nYear-1)], c(1, 2), mean)
  sYF.mean <- rowMeans(sYF[, 1:(nYear-1)], na.rm = T)
  sSA.mean <- rowMeans(sSA[, 1:(nYear-1)], na.rm = T)
  sAD.mean <- apply(sAD[, , 1:(nYear-1)], c(1, 2), mean)
  
  # population sizes
  nYF.mean <- rowMeans(nYF[, 1:nYear], na.rm = T)
  nSA.mean <- rowMeans(nSA[, 1:nYear], na.rm = T)
  nAD.mean <- apply(nAD[, , 1:(nYear-1)], c(1, 2), mean)
  nTOT.mean <- rowMeans(nTOT[, 1:nYear], na.rm = T)
  pYF.mean <- rowMeans(pYF[, 1:nYear], na.rm = T)
  pSA.mean <- rowMeans(pSA[, 1:nYear], na.rm = T)
  pAD.mean <- apply(pAD[, , 1:(nYear-1)], c(1, 2), mean)
  lambda.mean <- rowMeans(lambda[, 1:(nYear-1)], na.rm = T)
  
  
  ## Assemble list -------------------------------------------------------------
  
  paramSamples <- list(
    
    t = list(Bt = Bt,
             sPY = sPY,
             sYF = sYF,
             sSA = sSA,
             sAD = sAD,
             nYF = nYF,
             nSA = nSA,
             nAD = nAD,
             nTOT = nTOT,
             pYF = pYF,
             pSA = pSA,
             pAD = pAD,
             lambda = lambda),
    
    t.mean = list(Bt.mean = Bt.mean,
                  sPY.mean = sPY.mean,
                  sYF.mean = sYF.mean,
                  sSA.mean = sSA.mean,
                  sAD.mean = sAD.mean,
                  nYF.mean = nYF.mean,
                  nSA.mean = nSA.mean,
                  nAD.mean = nAD.mean,
                  nTOT.mean = nTOT.mean,
                  pYF.mean = pYF.mean,
                  pSA.mean = pSA.mean,
                  pAD.mean = pAD.mean,
                  lambda.mean = lambda.mean)
  )
  
  if(saveList){
    saveRDS(paramSamples, "results/paramSamples.rds")
  }
  
  return(paramSamples)
  
}

