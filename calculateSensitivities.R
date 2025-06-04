#' Calculate transient sensitivities & elasticities
#'
#' @param paramSamples list. Contains lists of posterior samples for all vital rates & population-level quantities.
#' @param nAge integer. Maximum age to consider in the analysis. nAge = 19 by default.
#' @param t.period vector. Optional argument specifying the years to use to calculate sensitivities & elasticities.
#'
#' @returns a list of lists containing posterior samples of transient sensitivities & elasticities for all vital rates & population structure (n).
#' @export
#'
#' @examples

calculateSensitivities <- function(paramSamples, nAge = 19, t.period = NULL){
  
  # # for testing purposes
  # # source('extractParamSamples.R')
  # # out.mcmc <- readRDS('results/IPM_CJSen_RSen_AB.rds')
  # # paramSamples <- extractParamSamples(MCMCsamples = out.mcmc, saveList = TRUE)
  # paramSamples <- readRDS('results/paramSamples.rds')
  # t.period <- NULL
  # nAge <- 19
  
  
  ## Calculate transient sensitivities -----------------------------------------
  
  # unpack objects from parameter list
  if(is.null(t.period)){
    for(i in 1:length(paramSamples$t.mean)){
      assign(names(paramSamples$t)[i], paramSamples$t.mean[[i]])
    }
  }else{
    for(i in 1:length(paramSamples$t.mean)){
      paramName <- names(paramSamples$t)[i]
      t.offset <- ifelse(paramName %in% c("Bt", "Ra", "sYAF", "sSA", "sAD"), 1, 0)
      
      focalParam <- paramSamples$t[[i]]
      
      if(length(dim(focalParam)) > 2){
        focalParam.mean <- apply(focalParam[, , t.period + t.offset], c(1, 2), mean)
      }else{
        focalParam.mean <- rowMeans(focalParam[, t.period + t.offset])
      }
      assign(paramName, focalParam.mean)
    }
  }
  
  # set sample number
  nSamples <- length(lambda)
  
  # set up list of arrays for storing transient sensitivities
  sensList <- list(
    sens.Bt = rep(NA, nSamples),
    sens.Ra = matrix(NA, nrow = nSamples, ncol = nAge),
    sens.sYAF = rep(NA, nSamples),
    sens.sSA = rep(NA, nSamples),
    sens.sAD = matrix(NA, nrow = nSamples, ncol = nAge),
    sens.pYAF = rep(NA, nSamples),
    sens.pSA = rep(NA, nSamples),
    sens.pAD = matrix(NA, nrow = nSamples, ncol = nAge)
  )
  
  # calculate transient sensitivities 
  # for vital rates & population size/structure (at the temporal mean)
  for(i in 1:nSamples){
    sensList$sens.Bt[i] <- sum(pAD[i, 2:nAge] * sAD[i, 2:nAge] * 0.5 * Ra[i, 2:nAge])
    
    sensList$sens.sYAF[i] <- pYAF[i]
    sensList$sens.sSA[i]  <- pSA[i]
    sensList$sens.pYAF[i] <- sYAF[i]
    sensList$sens.pSA[i]  <- sSA[i]
    
    for(a in 1:nAge){
      if(a == 1){
        sensList$sens.Ra[i, a] <- 0
        sensList$sens.sAD[i, a] <- 0
        sensList$sens.pAD[i, a] <- 0
      }else{
        sensList$sens.Ra[i, a] <- pAD[i, a] * sAD[i, a] * 0.5 * Bt[i]
        sensList$sens.sAD[i, a] <- pAD[i, a] * (1 + 0.5 * Bt[i] * Ra[i, a])
        sensList$sens.pAD[i, a] <- sAD[i, a] * (1 + 0.5 * Bt[i] * Ra[i, a])
      }
    }
  }
  
  # get posterior summaries for transient sensitivities
  postSum.sens <- data.frame()
  
  for(i in 1:length(sensList)){
    if(is.matrix(sensList[[i]])){
      
      # extract quantiles for age-specific parameters
      quantiles <- apply(sensList[[i]], 2, stats::quantile, probs = c(0.025, 0.5, 0.975))
      dimnames(quantiles)[[1]] <- c("lCI", "median", "uCI")
      data.temp <- cbind(data.frame(Parameter = names(sensList[i]), Age = 1:nAge), t(quantiles))
      
      # extract quantiles for sensitivities summed over ages
      sum.temp <- cbind(data.frame(Parameter = names(sensList[i]), Age = "summed"),
                        t(quantile(rowSums(sensList[[i]]), probs = c(0.025, 0.5, 0.975))))
      colnames(sum.temp)[3:5] <- c("lCI", "median", "uCI")
      data.temp <- rbind(data.temp, sum.temp)
      
      # merge into storage dataframe
      postSum.sens <- rbind(postSum.sens, data.temp)
    }else{
      
      # extract quantiles for age-specific parameters
      quantiles <- quantile(sensList[[i]], probs = c(0.025, 0.5, 0.975))
      names(quantiles) <- c("lCI", "median", "uCI")
      data.temp <- cbind(data.frame(Parameter = names(sensList[i]), Age = NA), t(quantiles))
      
      # merge into storage dataframe
      postSum.sens <- rbind(postSum.sens, data.temp)
    }
  }
  
  
  ## Calculate transient elasticities ------------------------------------------
  
  # calculate transient elasticities for vital rates & population size/structure
  # (evaluated at the temporal mean)
  elasList <- list(
    elas.Bt = sensList$sens.Bt * (Bt/lambda),
    # elas.Ra = sensList$sens.Ra * (Ra/lambda)
    elas.Ra = sensList$sens.Ra /lambda, # TODO: DISCUSS
    elas.sYAF = sensList$sens.sYAF * (sYAF/lambda),
    elas.sSA = sensList$sens.sSA * (sSA/lambda),
    elas.sAD = sensList$sens.sAD * (sAD/lambda),
    elas.pYAF = sensList$sens.pYAF * (pYAF/lambda),
    elas.pSA = sensList$sens.pSA * (pSA/lambda),
    elas.pAD = sensList$sens.pAD * (pAD/lambda)
  )
  
  # get posterior summaries for transient sensitivities
  postSum.elas <- data.frame()
  
  for(i in 1:length(elasList)){
    if(is.matrix(elasList[[i]])){
      
      # extract quantiles for age-specific parameters
      quantiles <- apply(elasList[[i]], 2, stats::quantile, probs = c(0.025, 0.5, 0.975))
      dimnames(quantiles)[[1]] <- c("lCI", "median", "uCI")
      data.temp <- cbind(data.frame(Parameter = names(elasList[i]), Age = 1:nAge), t(quantiles))
      
      # extract quantiles for elasticities summed over ages
      sum.temp <- cbind(data.frame(Parameter = names(elasList[i]), Age = "summed"),
                        t(quantile(rowSums(elasList[[i]]), probs = c(0.025, 0.5, 0.975))))
      colnames(sum.temp)[3:5] <- c("lCI", "median", "uCI")
      data.temp <- rbind(data.temp, sum.temp)
      
      # merge into storage dataframe
      postSum.elas <- rbind(postSum.elas, data.temp)
    }else{
      
      # extract quantiles for age-specific parameters
      quantiles <- quantile(elasList[[i]], probs = c(0.025, 0.5, 0.975))
      names(quantiles) <- c("lCI", "median", "uCI")
      data.temp <- cbind(data.frame(Parameter = names(elasList[i]), Age = NA), t(quantiles))
      
      # merge into storage dataframe
      postSum.elas <- rbind(postSum.elas, data.temp)
    }
  }
  
  
  ## Save & return results -----------------------------------------------------
  
  sensResults <- list(sensitivity = list(samples = sensList, summaries = postSum.sens),
                      elasticity = list(samples = elasList, summaries = postSum.elas))
  
  saveRDS(sensResults, file = "results/sensitivities.rds")
  return(sensResults)
  
}

