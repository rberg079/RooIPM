#' Run random design transient life table response experiment (LTRE)
#'
#' @param paramSamples list. Contains lists of posterior samples for all vital rates & population-level quantities.
#' @param sensitivities list. Contains lists containing posterior samples of transient sensitivities & elasticities for all vital rates & population structure (n).
#' @param nAge integer. Maximum age to consider in the analysis. nAge = 19 by default.
#' @param nYear integer. Number of time steps in the model. nYear = 17 by default.
#'
#' @returns a list of lists containing results of the LTRE analysis.
#' Object 'contList' is a list containing posterior distributions of all parameters' LTRE contributions (sublist 'cont'), as well as some auxiliary quantities (sublist 'other').
#' Object 'contData' is a dataframe containing posterior distributions for all parameters' LTRE contributions.
#' Object 'contData_summary' contains posterior summaries (medians & 95% credible intervals).
#' 
#' @export
#'
#' @examples

runLTRE_randomDesign <- function(paramSamples, sensitivities, nAge = 19, nYear = 17){
  
  # # for testing purposes
  # source('extractParamSamples.R')
  # source('calculateSensitivities.R')
  # out.mcmc <- readRDS('results/IPM_CJSen_RSen_AB.rds')
  # paramSamples <- extractParamSamples(MCMCsamples = out.mcmc, saveList = TRUE)
  # sensitivities <- calculateSensitivities(paramSamples = paramSamples)
  # 
  # OR
  # paramSamples <- readRDS('results/paramSamples.rds')
  # sensitivities <- readRDS('results/sensitivities.rds')
  # nYear <- 17
  # nAge <- 19
  
  
  ## Set up --------------------------------------------------------------------
  
  # load packages
  suppressPackageStartupMessages(library(tidyverse))
  
  # set sample number
  nSamples <- length(paramSamples$t.mean$lambda)
  
  dropParams <- c("nTOT")
  dropIdx <- which(names(paramSamples$t) %in% dropParams)
  
  paramList <- paramSamples$t[-dropIdx]
  sensList <- sensitivities$sensitivity$samples
  
  # set up list of arrays for storing calculated LTRE contributions
  contList <- list()
  
  for(x in 1:length(paramList)){
    if(names(paramList[x]) == "lambda"){
      next
    }
    if(length(dim(paramList[[x]])) == 3){
      tempList <- list()
      for(a in 1:nAge){
        tempList <- c(tempList, list(as.numeric(rep(NA, nSamples))))
      }
      names(tempList) <- paste0(names(paramList)[x], "_", 1:nAge)
    }else{
      tempList <- list(as.numeric(rep(NA, nSamples)))
      names(tempList) <- names(paramList)[x]
    }
    contList <- c(contList, tempList)
  }
  
  contCount <- length(contList)
  
  contList$est.var <- matrix(as.numeric(NA), nrow = contCount, ncol = nSamples)
  contList$est.covar <- matrix(as.numeric(NA), nrow = contCount, ncol = nSamples)
  
  
  ## Calculate LTRE contributions per sample (random design) -------------------
  
  for(i in 1:nSamples){
    
    # make lists of vital rates/population structure & of sensitivities
    dp.stoch.list <- list()
    sens.list <- list()
    
    for(x in 1:length(paramList)){
      if(names(paramList)[x] == "lambda"){
        next
      }
      
      # set time interval based on parameter
      if(names(paramList)[x] %in% c("nYAF", "nSA", "nAD")){
        tInt <- 2:nYear
      }else{
        tInt <- 1:(nYear-1)
      }
      
      # expand age if required & list relevant parameter estimates
      if(length(dim(paramList[[x]])) == 3){
        
        tempList <- list()
        tempListS <- list()
        
        for(a in 1:nAge){
          tempList <- c(tempList, list(paramList[[x]][i, a, tInt]))
          tempListS <- c(tempListS, list(sensList[[x]][i, a]))
        }
        names(tempList) <- paste0(names(paramList)[x], "_", 1:nAge)
        names(tempListS) <- paste0(names(paramList)[x], "_", 1:nAge)
      }else{
        tempList <- list(paramList[[x]][i, tInt])
        tempListS <- list(sensList[[x]][i])
        names(tempList) <- names(paramList)[x]
        names(tempListS) <- names(paramList)[x]
      }
      dp.stoch.list <- c(dp.stoch.list, tempList)
      sens.list <- c(sens.list, tempListS)
    }
    
    # convert parameter list to matrix
    dp.stoch <- as.matrix(dplyr::bind_rows(dp.stoch.list, .id = "column.label"))
    
    # derive process variances & covariances
    dp.varcov <- var(dp.stoch)
    
    # save total estimated (co)variance per parameter
    contList$est.var[,i] <- diag(dp.varcov)
    contList$est.covar[,i] <- rowSums(dp.varcov, na.rm = T)
    
    # convert sensitivity list to vector
    sens.vec <- do.call(c, sens.list)
    
    # calculate demographic contributions
    # NOTE: here we multiply sensitivities & (co)variances
    cont.mat <- matrix(NA, nrow = length(sens.vec), ncol = length(sens.vec))
    for(k in 1:length(sens.vec)){
      for(l in 1:length(sens.vec)){
        cont.mat[k, l] <- dp.varcov[k, l] * sens.vec[k] * sens.vec[l]
      }
    }
    
    # summarise contributions (sum of variances & covariances)
    cont <- rowSums(cont.mat)
    names(cont) <- names(sens.vec)
    
    # insert contributions into storage list
    for(x in 1:length(cont)){
      contList[[x]][i] <- cont[x]
    }
  }
  
  # restructure results list
  contList <- list(cont = contList[1:contCount],
                   other = list(est.var = contList$est.var,
                                est.covar = contList$est.covar)
  )
  
  
  ## Summarise results ---------------------------------------------------------
  
  # check sum of contributions against variance in lambda
  contList$other$total.contSum <- rowSums(dplyr::bind_rows(contList$cont))
  quantile(contList$other$total.contSum, probs = c(0.025, 0.5, 0.975))
  
  contList$other$tempvar.lambda <- matrixStats::rowVars(paramList$lambda[, 1:(nYear-1)])
  quantile(contList$other$tempvar.lambda, probs = c(0.025, 0.5, 0.975))
  
  # calculate summed contributions for age-specific parameters
  sumParams <- names(paramList)[which(names(paramList) %in% c("nAD", "sAD", "Ra"))]
  
  for(x in 1:length(sumParams)){
    subList <- contList$cont[which(grepl(paste0(sumParams[x], "_"), names(contList$cont)))]
    
    contList$cont$newSum <- rowSums(dplyr::bind_rows(subList))
    names(contList$cont)[which(names(contList$cont) == "newSum")] <- paste0(sumParams[x], "_sum")
  }
  
  # arrange results as dataframe
  # contData <- melt(dplyr::bind_rows(contList$cont, .id = "column.label"))
  contData <- contList$cont %>% 
    dplyr::bind_rows() %>% 
    tidyr::pivot_longer(cols = everything())
  
  contData <- cbind(contData, stringr::str_split_fixed(contData$name, "_", 2))
  colnames(contData) <- c("Variable", "Contribution", "Parameter", "Age")
  
  contData$Age[which(contData$Age == "")] <- NA
  
  # make posterior summaries
  contData_summary <- contData %>% 
    dplyr::group_by(Variable) %>% 
    dplyr::summarise(lCI = quantile(Contribution, 0.025),
                     median = median(Contribution),
                     uCI = quantile(Contribution, 0.975))
  
  
  ## Collect & return results --------------------------------------------------
  
  results <- list(contList = contList,
                  contData = contData,
                  contData_summary = contData_summary)
  
  saveRDS(results, file = "results/LTREresults_random.rds")
  
  return(results)
  
}

