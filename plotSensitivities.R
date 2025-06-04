#' Plot transient sensitivities & elasticities
#'
#' @param sensitivities list. Contains lists containing posterior samples of transient sensitivities & elasticities for all vital rates & population structure (n).
#' @param nAge integer. Maximum age to consider in the analysis. nAge = 19 by default.
#' @param plotFolder character string. Path to the folder in which to store plots.
#'
#' @returns a character vector of plot names. The plots themselves are saved as pdfs in plotFolder.
#' @export
#'
#' @examples

plotSensitivities <- function(sensitivities, nAge = 19, plotFolder){

  # for testing purposes
  sensitivities <- readRDS('results/sensitivities.rds')
  plotFolder = c("figures")
  nAge = 19
  
  
  ## Set up --------------------------------------------------------------------
  
  library(coda)
  suppressPackageStartupMessages(library(tidyverse))
  suppressPackageStartupMessages(library(data.table))
  # library(NatParksPalettes)
  library(paletteer)
  
  # make plotting directory if it does not exist already
  if(!dir.exists(plotFolder)){
    dir.create(plotFolder)
  }
  
  for(i in 1:2){
    
    # select relevant data
    params <- sensitivities[[i]]$samples
    
    # drop pre-fix for generalising
    names(params) <- sub("^.{4}\\.", "", names(params))
    
    # extract number of samples
    nSamples <- length(params[[1]])
    
    ## Assemble summarised data ------------------------------------------------
    
    # assemble data
    sum.data <- data.frame(
      type = rep(c("Birth rate",
                   "Survival of jellybeans",
                   "Survival of young-at-foot",
                   "Survival of subadults",
                   "Survival of adults",
                   "Proportion of young-at-foot",
                   "Proportion of subadults",
                   "Proportion of adults")),
      estimate = c(params$Bt,
                   rowSums(params$Ra),
                   params$sYAF,
                   params$sSA,
                   rowSums(params$sAD),
                   params$pYAF,
                   params$pSA,
                   rowSums(params$pAD)))
    
    # order factor levels
    sum.data$type <- factor(sum.data$type,
                            levels = c("Birth rate",
                                       "Survival of jellybeans",
                                       "Survival of young-at-foot",
                                       "Survival of subadults",
                                       "Survival of adults",
                                       "Proportion of young-at-foot",
                                       "Proportion of subadults",
                                       "Proportion of adults"))
    
    
    ## Assemble age-specific data ----------------------------------------------
    
    # bind all data into a dataframe
    age.data <- data.frame(rlist::list.cbind(params))
    
    # change column names
    colnames(age.data) <- c("Bt", paste0("Ra_", 1:nAge),
                            "sYAF", "sSA", paste0("sAD_", 1:nAge),
                            "pYAF", "pSA", paste0("pAD_", 1:nAge))
    
    # convert to longitudinal format
    age.data <- reshape2::melt(age.data)
    
    
    ## Plot sensitivities/elasticities -----------------------------------------
    
    # plot colours
    temp.colours <- paletteer::paletteer_c("grDevices::Temps", length(unique(sum.data$type))) # -3?
    # plot.colours <- c("#047993FF", "#005F94FF", temp.colours[1:3], rep(temp.colours[4], 2), temp.colours[5:6])
    plot.colours <- temp.colours
    
    # summed estimates for all parameters
    addline_format <- function(x,...){
      gsub('\\s','\n',x)
    }
    
    p.sum <- ggplot(sum.data, aes(x = type, y = estimate, group = type)) + 
      geom_violin(aes(fill = type, colour = type), alpha = 0.5, scale = 'width', draw_quantiles = 0.5) + 
      geom_hline(yintercept = 0, colour = "grey70", linetype = "dashed") + 
      ylab(ifelse(i == 1, "Sensitivity", "Elasticity")) + 
      xlab('') + 
      scale_fill_manual(values = plot.colours) + 
      scale_colour_manual(values = plot.colours) + 
      scale_x_discrete(labels = c("Birth\nrate",
                                  "Surv. of\nbeans",
                                  "Surv. of\nYAFs",
                                  "Surv. of\nsubadults",
                                  "Surv. of\nadults",
                                  "Prop. of\nYAFs",
                                  "Prop. of\nsubadults",
                                  "Prop. of\nadults")) + 
      theme_bw() + 
      theme(legend.position = 'none',
            panel.grid = element_blank(),
            axis.text.x = element_text(size = 10),
            axis.title = element_text(size = 10))
    
    # reproductive success panel
    R.colours <- c(plot.colours[1], rep(plot.colours[2], 18))
    names(R.colours) <- c("Bt", paste0("Ra_", 2:nAge))
    
    p.R <- ggplot(subset(age.data, variable %in% c("Bt", paste0("Ra_", 2:nAge)))) +
      geom_violin(aes(x = factor(variable, levels = c("Bt", paste0("Ra_", 2:nAge))),
                      y = value, fill = variable), alpha = 0.5, scale = 'width', draw_quantiles = 0.5, position = 'dodge') +
      ylab(ifelse(i == 1, "Sensitivity", "Elasticity")) + 
      xlab('') + 
      scale_x_discrete(labels = expression(B, R[2], R[3], R[4], R[5], R[6], R[7], R[8], R[9], R[10],
                                           R[11], R[12], R[13], R[14], R[15], R[16], R[17], R[18], R[19])) +
      scale_fill_manual(values = R.colours) +
      theme_bw() + 
      theme(legend.position = 'none',
            panel.grid = element_blank(), 
            axis.text.x = element_text(size = 10), 
            axis.title = element_text(size = 10))
    
    # survival panel
    S.colours <- c(plot.colours[3:4], rep(plot.colours[5], 18))
    names(S.colours) <- c("sYAF", "sSA", paste0("sAD_", 2:nAge))
    
    p.S <- ggplot(subset(age.data, variable %in% c("sYAF", "sSA", paste0("sAD_", 2:nAge)))) +
      geom_violin(aes(x = factor(variable, levels = c("sYAF", "sSA", paste0("sAD_", 2:nAge))),
                      y = value, fill = variable), alpha = 0.5, scale = 'width', draw_quantiles = 0.5, position = 'dodge') + 
      ylab(ifelse(i == 1, "Sensitivity", "Elasticity")) + 
      xlab('') + 
      scale_x_discrete(labels = expression(S[0], S[1],
                                           S[2], S[3], S[4], S[5], S[6], 
                                           S[7], S[8], S[9], S[10], S[11], S[12],
                                           S[13], S[14], S[15], S[16], S[17], S[18], S[19])) +
      scale_fill_manual(values = S.colours) +
      theme_bw() + 
      theme(legend.position = 'none',
            panel.grid = element_blank(), 
            axis.text.x = element_text(size = 10), 
            axis.title = element_text(size = 10))
    
    # population structure panel
    P.colours <- c(plot.colours[6:7], rep(plot.colours[8], 18))
    names(P.colours) <- c("pYAF", "pSA", paste0("pAD_", 2:nAge))
    
    p.P <- ggplot(subset(age.data, variable %in% c("pYAF", "pSA", paste0("pAD_", 2:nAge)))) +
      geom_violin(aes(x = factor(variable, levels = c("pYAF", "pSA", paste0("pAD_", 2:nAge))),
                      y = value, fill = variable), alpha = 0.5, scale = 'width', draw_quantiles = 0.5, position = 'dodge') +
      ylab(ifelse(i == 1, "Sensitivity", "Elasticity")) + 
      xlab('') + 
      scale_x_discrete(labels = expression(P[0], P[1],
                                           P[2], P[3], P[4], P[5], P[6], 
                                           P[7], P[8], P[9], P[10], P[11], P[12],
                                           P[13], P[14], P[15], P[16], P[17], P[18], P[19])) +
      scale_fill_manual(values = P.colours) +
      theme_bw() + 
      theme(legend.position = 'none',
            panel.grid = element_blank(), 
            axis.text.x = element_text(size = 10), 
            axis.title = element_text(size = 10))
    
    # combine panels & save to pdf
    pdf(paste0(plotFolder, ifelse(i == 1, "/sensitivities", "/elasticities"), "_sum.pdf"), width = 8, height = 4) # or 10 & 6
    print(
      p.sum
    )
    dev.off()
    
    library(patchwork)
    pdf(paste0(plotFolder, ifelse(i == 1, "/sensitivities", "/elasticities"), "_age.pdf"), width = 12, height = 6) # or 10 & 6
    print(
      # p.sum + (p.R / p.S / p.P)
      (p.sum + labs(tag = 'a)')) + ((p.R + labs(tag = 'b)')) / (p.S + labs(tag = 'c)')) / (p.P + labs(tag = 'd)')))
    )
    dev.off()
  }
  
  
  ## Return list of plots ------------------------------------------------------
  
  plotList <- c(paste0(plotFolder, c("Sensitivities", "Elasticities"), "_sum.pdf"),
                paste0(plotFolder, c("Sensitivities", "Elasticities"), "_age.pdf"))
  
  return(plotList)
  
}

