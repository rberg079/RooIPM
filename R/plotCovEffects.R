#' Plot vital rates as a function of covariates
#'
#' @param MCMCsamples mcmc.list or matrix containing the output of the fitted IPM. 
#' @param ages integer vector. Age classes to plot. ages = c(2, 6, 10, 14) by default.
#' @param densMIN numeric. Minimum value of standardized density covariate. densMIN = -2 by default.
#' @param densMAX numeric. Maximum value of standardized density covariate. densMAX = 2 by default.
#' @param vegMIN numeric. Minimum value of standardized vegetation covariate. vegMIN = -2 by default.
#' @param vegMAX numeric. Maximum value of standardized vegetation covariate. vegMAx = 2 by default.
#' @param plotFolder character string. Path to the folder in which to store plots.
#'
#' @return character vector of plot names.
#' The plots themselves are saved as pdf's in the specified subfolder.
#' @export
#'
#' @examples

plotIPM_DemographicCovariates <- function(MCMCsamples, 
                                          ages = c(1, 4, 7, 10), 
                                          densMIN = -2, densMAX = 2, 
                                          vegMIN = -2, vegMAX = 2,
                                          plotFolder){
  
  # # for testing purposes
  # MCMCsamples = readRDS('results/IPM_CJSen_RSen_AB_DynDens_dCJS_12_noW_25BR_S33_R22_stochV_8chains.rds')
  # ages = c(1, 4, 7, 10)
  # densMIN = -2
  # densMAX = 2
  # vegMIN = -2
  # vegMAX = 2
  # plotFolder = "figures/densityChecks/varyNvarsSR"
  
  
  ## Set up --------------------------------------------------------------------
  
  suppressPackageStartupMessages(library(tidyverse))
  library(coda)
  
  # create plot folder if not present
  if(!dir.exists(plotFolder)){dir.create(plotFolder)}
  
  # convert MCMC samples to matrix
  out.mat <- as.matrix(MCMCsamples)
  
  # make vectors of covariate values for predictions
  densCov <- seq(densMIN, densMAX, length.out = 100)
  vegCov  <- seq(vegMIN, vegMAX, length.out = 100)
  
  # set up empty dataframe to store results
  pred.data <- data.frame()
  
  # helper function to check for age-specific slopes vs general slopes
  get_slope <- function(out_matrix, specific_param, general_param){
    if(specific_param %in% colnames(out_matrix)){
      return(out_matrix[, specific_param])
    }else if(general_param %in% colnames(out_matrix)){
      return(out_matrix[, general_param])
    }else{
      return(rep(0, nrow(out_matrix)))
    }
  }
  
  
  ## Density predictions -------------------------------------------------------
  
  for(x in 1:length(densCov)){
    for(a in ages){
      
      # determine correct density slopes based on age
      if(a == 1){
        bD.S <- get_slope(out.mat, "BetaD.Sy", "BetaD.S")
        bD.R <- rep(0, nrow(out.mat)) # young don't reproduce
      }else if(a >= 2 & a <= 9){
        bD.S <- get_slope(out.mat, "BetaD.Sp", "BetaD.S")
        bD.R <- get_slope(out.mat, "BetaD.Rp", "BetaD.R")
      }else{
        bD.S <- get_slope(out.mat, "BetaD.So", "BetaD.S")
        bD.R <- get_slope(out.mat, "BetaD.Ro", "BetaD.R")
      }
      
      # survival prediction
      S.pred <- plogis(qlogis(out.mat[, paste0("Mu.S[", a, "]")]) + 
                         bD.S * densCov[x])
      
      # survival of pouch young prediction
      if(a == 1){
        sPY.pred <- rep(0, nrow(out.mat))
      }else{
        sPY.pred <- plogis(qlogis(out.mat[, paste0("Mu.R[", a, "]")]) + 
                             bD.R * densCov[x])
      }
      
      # assemble in temporary data frame
      temp <- rbind(quantile(S.pred, probs = c(0.025, 0.5, 0.975)),
                    quantile(sPY.pred, probs = c(0.025, 0.5, 0.975))) %>%
        as.data.frame() %>%
        dplyr::rename(lCI = 1, median = 2, uCI = 3) %>%
        dplyr::mutate(VitalRate = c("Survival (S)", "Survival of pouch young (sPY)"),
                      Age = as.factor(a),
                      CovariateValue = densCov[x],
                      CovariateType = "Density")
      
      # combine
      pred.data <- rbind(pred.data, temp)
    }
  }
  
  
  ## Vegetation predictions ----------------------------------------------------
  
  for(x in 1:length(vegCov)){
    for(a in ages){
      
      # determine correct vegetation slopes based on age
      if(a == 1){
        bV.S <- get_slope(out.mat, "BetaV.Sy", "BetaV.S")
      }else if(a >= 2 & a <= 9){
        bV.S <- get_slope(out.mat, "BetaV.Sp", "BetaV.S")
      }else{
        bV.S <- get_slope(out.mat, "BetaV.So", "BetaV.S")
      }
      
      # survival prediction
      S.pred <- plogis(qlogis(out.mat[, paste0("Mu.S[", a, "]")]) + 
                         bV.S * vegCov[x])
      
      # assemble in temporary data frame (sPY removed for vegetation)
      temp <- quantile(S.pred, probs = c(0.025, 0.5, 0.975)) %>%
        t() %>%
        as.data.frame() %>%
        dplyr::rename(lCI = 1, median = 2, uCI = 3) %>%
        dplyr::mutate(VitalRate = "Survival (S)",
                      Age = as.factor(a),
                      CovariateValue = vegCov[x],
                      CovariateType = "Vegetation")
      
      # combine
      pred.data <- rbind(pred.data, temp)
    }
  }
  
  
  ## Plot predictions ----------------------------------------------------------
  
  # re-arrange factor levels
  pred.data$VitalRate <- factor(pred.data$VitalRate, levels = c("Survival (S)", "Survival of pouch young (sPY)"))
  
  # dynamically assign colors based on the number of ages provided
  default_cols <- c("#B3B3B3", "#52BA88FF", "#089392FF", "#BED68AFF", "#E79069FF", "#D9A066", "#8F6449")
  age.cols <- setNames(default_cols[1:length(ages)], as.character(ages))
  
  # density effects plot (4x2 grid via facet_grid)
  p.dens <- pred.data %>%
    dplyr::filter(CovariateType == "Density") %>%
    ggplot(aes(x = CovariateValue, group = Age)) + 
    geom_line(aes(y = median, color = Age)) + 
    geom_ribbon(aes(ymin = lCI, ymax = uCI, fill = Age), alpha = 0.3) + 
    scale_color_manual(values = age.cols) + 
    scale_fill_manual(values = age.cols) + 
    ylab("Predicted probability") + 
    xlab("Population density (standardized)") + 
    facet_grid(VitalRate ~ Age, scales = "free_y") + 
    theme_bw() + 
    theme(panel.grid = element_blank(), legend.position = "none")
  
  pdf(paste0(plotFolder, "/densityEffects.pdf"), width = 8, height = 5)
  print(p.dens)
  dev.off()
  
  # vegetation effects plot (4x1 grid via facet_grid, since sPY is excluded)
  p.veg <- pred.data %>%
    dplyr::filter(CovariateType == "Vegetation") %>%
    ggplot(aes(x = CovariateValue, group = Age)) + 
    geom_line(aes(y = median, color = Age)) + 
    geom_ribbon(aes(ymin = lCI, ymax = uCI, fill = Age), alpha = 0.3) + 
    scale_color_manual(values = age.cols) + 
    scale_fill_manual(values = age.cols) + 
    ylab("Predicted probability") + 
    xlab("Vegetation (standardized)") + 
    facet_grid(VitalRate ~ Age, scales = "free_y") + 
    theme_bw() + 
    theme(panel.grid = element_blank(), legend.position = "none")
  
  pdf(paste0(plotFolder, "/vegetationEffects.pdf"), width = 8, height = 3)
  print(p.veg)
  dev.off()
  
  # return list of plots
  plotList <- c(paste0(plotFolder, "/densityEffects.pdf"),
                paste0(plotFolder, "/vegetationEffects.pdf"))
  
  return(plotList)
}

