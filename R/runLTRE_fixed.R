#' Run fixed design transient life table response experiment (LTRE) for a pair of years
#'
#' @param paramSamples list. Contains lists of posterior samples for all vital rates & population-level quantities.
#' @param t.pair vector. Contains 2 integers specifying the indices of the two years to compare in the analysis.
#' @param nAge integer. Maximum age to consider in the analysis. nAge = 19 by default.
#'
#' @returns a list of lists containing results of the LTRE analysis.
#' Object 'contList' is a list containing posterior distributions of all parameters' LTRE contributions (sublist 'cont'), as well as some auxiliary quantities (sublist 'other').
#' Object 'contData' is a dataframe containing posterior distributions for all parameters' LTRE contributions.
#' Object 'contData_summary' contains posterior summaries (medians & 95% credible intervals).
#' @export
#'
#' @examples

runLTRE_fixed <- function(paramSamples, t.pair, nAge = 19){
  
  # # for testing purposes
  # source('extractParamSamples.R')
  # out.mcmc <- readRDS('results/IPM_CJSen_RSen_AB.rds')
  # paramSamples <- extractParamSamples(MCMCsamples = out.mcmc, saveList = TRUE)
  # 
  # OR
  paramSamples <- readRDS('results/paramSamples.rds')
  t.pair <- c(4, 8)
  nAge <- 19
  
  
  ## Set up --------------------------------------------------------------------
  
  # load packages
  suppressPackageStartupMessages(library(tidyverse))
  
  # set sample number
  nSamples <- length(paramSamples$t.mean$lambda)
  
  dropParams <- c("nTOT")
  dropIdx <- which(names(paramSamples$t) %in% dropParams)
  
  paramList <- paramSamples$t[-dropIdx]
  
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
  
  contList$delta_lambda <- rep(NA, nSamples)
  contList$delta_loglambda <- rep(NA, nSamples)
  
  
  ## Calculate sensitivities for relevant year pair ----------------------------
  
  source('calculateSensitivities.R')
  sensitivities <- calculateSensitivities(paramSamples = paramSamples,
                                          nAge = nAge,
                                          t.period = t.pair)
  
  sensList <- sensitivities$sensitivity$samples
  
  
  ## Calculate LTRE contributions per sample (fixed design) --------------------
  
  for(i in 1:nSamples){
    
    # make lists of vital rates/population structure & of sensitivities
    param_list <- list()
    sens_list <- list()
    
    for(x in 1:length(paramList)){
      if(names(paramList)[x] == "lambda"){
        next
      }
      
      # set time interval based on parameter
      if(names(paramList)[x] %in% c("nYAF", "nSA", "nAD")){
        tInt <- t.pair + 1
      }else{
        tInt <- t.pair
      }
      
      # expand age if required & list relevant parameter estimates
      if(length(dim(paramList[[x]])) == 3){
        
        tempList <- list()
        tempListS <- list()
        
        for(a in 1:nAge){
          tempList <- c(tempList, list(diff(paramList[[x]][i, a, tInt])))
          tempListS <- c(tempListS, list(sensList[[x]][i, a]))
        }
        names(tempList) <- paste0(names(paramList)[x], "_", 1:nAge)
        names(tempListS) <- paste0(names(paramList)[x], "_", 1:nAge)
      }else{
        tempList <- list(diff(paramList[[x]][i, tInt]))
        tempListS <- list(sensList[[x]][i])
        names(tempList) <- names(paramList)[x]
        names(tempListS) <- names(paramList)[x]
      }
      param_list <- c(param_list, tempList)
      sens_list <- c(sens_list, tempListS)
    }
    
    # convert parameter & sensitivity lists to vectors
    paramvec <- do.call(c, param_list)
    sensvec <- do.call(c, sens_list)
    
    # calculate & save change in lambda
    contList$delta_lambda[i] <- diff(paramList$lambda[i, t.pair])
    contList$delta_loglambda[i] <- diff(log(paramList$lambda[i, t.pair]))
    
    # calculate demographic contributions
    # NOTE: here we multiply sensitivities & parameter differences
    cont <- paramvec * sensvec
    names(cont) <- names(sensvec)
    
    # insert contributions into storage list
    for(x in 1:length(cont)){
      contList[[x]][i] <- cont[x]
    }
  }
  
  # restructure results list
  contList <- list(cont = contList[1:contCount],
                   other = list(delta_lambda = contList$delta_lambda,
                                delta_loglambda = contList$delta_loglambda)
  )
  
  
  ## Summarise results ---------------------------------------------------------
  
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
  
  saveRDS(results, file = "results/LTREresults_fixed.rds")
  
  return(results)
  
}

