#' Plots results from a random design transient LTRE
#'
#' @param LTREresults list. Contains results of a random design transient LTRE.
#' @param nAge integer. Maximum age to consider in the analysis. nAge = 19 by default.
#' @param plotFolder character string. Path to the folder in which to store plots.
#'
#' @returns a character vector of plot names. The plots themselves are saved as pdfs in plotFolder.
#' @export
#'
#' @examples

plotLTRE_random <- function(LTREresults, nAge = 19, plotFolder){
  
  # # for testing purposes
  # LTREresults <- readRDS('results/LTREresults_random.rds')
  # plotFolder = c("figures")
  # nAge = 19
  
  
  ## Set up --------------------------------------------------------------------
  
  # load packages
  suppressPackageStartupMessages(library(tidyverse))
  
  # make plotting directory if it does not exist already
  if(!dir.exists(plotFolder)){
    dir.create(plotFolder)
  }
  
  
  ## Format data ---------------------------------------------------------------
  
  # select relevant data
  contData <- LTREresults$contData
  
  # extract number of samples
  nSamples <- length(LTREresults$contList$cont[[1]])
  
  # split & format summed data
  contData_sum <- contData %>% 
    dplyr::filter(Variable %in% c("Bt", "Ra_sum",
                                  "sYAF", "sSA", "sAD_sum",
                                  "pYAF", "pSA", "pAD_sum")) %>%
    dplyr::mutate(type = dplyr::case_when(Variable == "Bt" ~ "Birth rate",
                                          Variable == "Ra_sum" ~ "Survival of jellybeans",
                                          Variable == "sYAF" ~ "Survival of young-at-foot",
                                          Variable == "sSA" ~ "Survival of subadults",
                                          Variable == "sAD_sum" ~ "Survival of adults",
                                          Variable == "pYAF" ~ "Proportion of young-at-foot",
                                          Variable == "pSA" ~ "Proportion of subadults",
                                          Variable == "pAD_sum" ~ "Proportion of adults"))
  
  # make ordered list of parameter types
  typeList <- c("Birth rate", "Survival of jellybeans",
                "Survival of young-at-foot", "Survival of subadults", "Survival of adults",
                "Proportion of young-at-foot", "Proportion of subadults", "Proportion of adults")
  
  # order factor levels
  contData_sum$type <- factor(contData_sum$type, levels = typeList)
  
  
  ## Plot contributions (violin plots) -----------------------------------------
  
  # plot colours
  temp.colours <- paletteer::paletteer_c("grDevices::Temps", length(unique(contData_sum$type)))
  # plot.colours <- c("#047993FF", "#005F94FF", temp.colours[1:3], rep(temp.colours[4], 2), temp.colours[5:6])
  plot.colours <- temp.colours
  
  # summed estimates for all parameters
  addline_format <- function(x,...){
    gsub('\\s','\n',x)
  }
  
  p.sum <- contData_sum %>% 
    # dplyr::filter(Contribution < 0.5) %>% 
    ggplot(aes(x = type, y = Contribution, group = type)) +
    geom_violin(aes(fill = type, colour = type), alpha = 0.5, scale = "width", draw_quantiles = 0.5) +
    geom_hline(yintercept = 0, colour = "grey70", linetype = "dashed") +
    ylab("Contribution") +
    xlab("") +
    scale_fill_manual(values = plot.colours) +
    scale_colour_manual(values = plot.colours) +
    # scale_x_discrete(labels = addline_format(typeList)) +
    scale_x_discrete(labels = c("Birth\nrate",
                                "Surv. of\nbeans",
                                "Surv. of\nYAFs",
                                "Surv. of\nsubadults",
                                "Surv. of\nadults",
                                "Prop. of\nYAFs",
                                "Prop. of\nsubadults",
                                "Prop. of\nadults")) +
    theme_bw() +
    theme(legend.position = "none",
          panel.grid = element_blank(),
          axis.text.x = element_text(size = 10),
          axis.title = element_text(size = 10))
  
  # reproductive success panel
  R.colours <- c(plot.colours[1], rep(plot.colours[2], 18))
  names(R.colours) <- c("Bt", paste0("Ra_", 2:nAge))
  
  p.R <- ggplot(subset(contData, Variable %in% c("Bt", paste0("Ra_", 2:nAge)))) +
    geom_violin(aes(x = factor(Variable, levels = c("Bt", paste0("Ra_", 2:nAge))),
                    y = Contribution, fill = Variable), alpha = 0.5, scale = "width", draw_quantiles = 0.5) +
    geom_hline(yintercept = 0, colour = "grey70", linetype = "dashed") +
    ylab("Contribution") +
    xlab("") +
    scale_x_discrete(labels = expression(B, R[2], R[3], R[4], R[5], R[6], R[7], R[8], R[9], R[10],
                                         R[11], R[12], R[13], R[14], R[15], R[16], R[17], R[18], R[19])) +
    scale_fill_manual(values = R.colours) +
    theme_bw() +
    theme(legend.position = "none",
          panel.grid = element_blank(),
          axis.text.x = element_text(size = 10),
          axis.title = element_text(size = 10))
  
  # survival panel
  S.colours <- c(plot.colours[3:4], rep(plot.colours[5], 18))
  names(S.colours) <- c("sYAF", "sSA", paste0("sAD_", 2:nAge))
  
  p.S <- ggplot(subset(contData, Variable %in% c("sYAF", "sSA", paste0("sAD_", 2:nAge)))) +
    geom_violin(aes(x = factor(Variable, levels = c("sYAF", "sSA", paste0("sAD_", 2:nAge))),
                    y = Contribution, fill = Variable), alpha = 0.5, scale = "width", draw_quantiles = 0.5) +
    geom_hline(yintercept = 0, colour = "grey70", linetype = "dashed") +
    ylab("Contribution") +
    xlab("") +
    scale_x_discrete(labels = expression(S[0], S[1],
                                         S[2], S[3], S[4], S[5], S[6], 
                                         S[7], S[8], S[9], S[10], S[11], S[12],
                                         S[13], S[14], S[15], S[16], S[17], S[18], S[19])) +
    scale_fill_manual(values = S.colours) +
    theme_bw() +
    theme(legend.position = "none",
          panel.grid = element_blank(),
          axis.text.x = element_text(size = 10),
          axis.title = element_text(size = 10))
  
  # population structure panel
  P.colours <- c(plot.colours[6:7], rep(plot.colours[8], 18))
  names(P.colours) <- c("nYAF", "nSA", paste0("nAD_", 2:nAge))
  
  p.P <- ggplot(subset(contData, Variable %in% c("nYAF", "nSA", paste0("nAD_", 2:nAge)))) +
    geom_violin(aes(x = factor(Variable, levels = c("nYAF", "nSA", paste0("nAD_", 2:nAge))),
                    y = Contribution, fill = Variable), alpha = 0.5, scale = "width", draw_quantiles = 0.5) +
    geom_hline(yintercept = 0, colour = "grey70", linetype = "dashed") +
    ylab("Contribution") +
    xlab("") +
    scale_x_discrete(labels = expression(P[0], P[1],
                                         P[2], P[3], P[4], P[5], P[6], 
                                         P[7], P[8], P[9], P[10], P[11], P[12],
                                         P[13], P[14], P[15], P[16], P[17], P[18], P[19])) +
    scale_fill_manual(values = P.colours) +
    theme_bw() +
    theme(legend.position = "none",
          panel.grid = element_blank(),
          axis.text.x = element_text(size = 10),
          axis.title = element_text(size = 10))
  
  # combine panels & save to pdf
  pdf(paste0(plotFolder, "/LTRErandom_sum.pdf"), width = 10, height = 6)
  print(
    p.sum
  )
  dev.off()
  
  library(patchwork)
  pdf(paste0(plotFolder, "/LTRErandom_age.pdf"), width = 7, height = 8)
  print(
    # p.sum + (p.R / p.S / p.P)
    (p.sum + labs(tag = "a)")) + ((p.R + labs(tag = "b)")) / (p.S + labs(tag = "c)")) / (p.P + labs(tag = "d)")))
  )
  dev.off()
  
  
  ## Return list of plots ------------------------------------------------------
  plotList <- c(paste0(plotFolder, "/LTRErandom_sum.pdf"),
                paste0(plotFolder, "/LTRErandom_age.pdf"))
  
  return(plotList)
  
}

