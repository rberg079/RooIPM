#' Plots results from a random design transient LTRE
#'
#' @param LTREresults list. Contains results of a random design transient LTRE.
#' @param nAge integer. Maximum age to consider in the analysis. nAge = 19 by default.
#' @param plotFolder character string. Path to the folder in which to store plots.
#' @param oneProp logical. If TRUE, sums age-structure proportions into a single age structure category. oneProp = TRUE by default.
#' @param returnSummary logical. If TRUE, returns saved file paths and summary data frame.. returnSummary = TRUE by default.
#'
#' @returns a character vector of plot names & optionally the LTRE summary metrics. The plots themselves are saved as pdfs in plotFolder.
#' @export
#'
#' @examples

plotLTRE_random <- function(LTREresults, nAge = 19, plotFolder,
                            oneProp = TRUE, returnSummary = TRUE){
  
  # # for testing purposes
  # LTREresults <- readRDS('results/LTREresults_random.rds')
  # plotFolder = c("figures/resultsDave2Covs")
  # nAge = 19
  # oneProp <- TRUE
  # returnSummary <- TRUE
  
  
  ## Set up --------------------------------------------------------------------
  
  # load packages
  suppressPackageStartupMessages(library(tidyverse))
  suppressPackageStartupMessages(library(patchwork))
  library(paletteer)
  library(scales)
  
  # make plotting directory if it does not exist already
  if(!dir.exists(plotFolder)){
    dir.create(plotFolder)
  }
  
  
  ## Assemble data -------------------------------------------------------------
  
  # select relevant data
  contData <- LTREresults$contData
  
  # assign draw indices to samples
  contData <- contData %>%
    mutate(draw = rep(1:(nrow(.) / length(unique(Variable))),
                      each = length(unique(Variable))))
  
  # split & format summed data
  if(oneProp) {
    
    tmp <- contData %>% 
      filter(Variable %in% c("pYF", "pSA", "pAD_all")) %>% 
      group_by(draw) %>% 
      summarise(Contribution = sum(Contribution), .groups = "drop") %>% 
      mutate(Variable = "p_all", type = "Age structure")
    
    plotData <- contData %>% 
      filter(Variable %in% c("BR_all", "sPY_all",
                             "sYF", "sSA", "sAD_2to9", "sAD_10up")) %>%
      mutate(type = case_when(
        Variable == "BR_all"   ~ "Birth rate",
        Variable == "sPY_all"  ~ "Survival of pouch young",
        Variable == "sYF"      ~ "Survival of young-at-foot",
        Variable == "sSA"      ~ "Survival of subadults",
        Variable == "sAD_2to9" ~ "Survival of adults (2–9)",
        Variable == "sAD_10up" ~ "Survival of adults (10+)"
      )) %>% 
      bind_rows(tmp)
    
    types <- c("Birth rate",
               "Survival of pouch young",
               "Survival of young-at-foot", 
               "Survival of subadults",
               "Survival of adults (2–9)",
               "Survival of adults (10+)",
               "Age structure")
    
    names <- c("Birth\nrate",
               "Survival of\npouch young",
               "Survival of\nyoung-at-foot", 
               "Survival of\nsubadults",
               "Survival of\nadults (2–9)",
               "Survival of\nadults (10+)",
               "Age\nstructure")
  }else{
    
    plotData <- contData %>% 
      filter(Variable %in% c("BR_all", "sPY_all",
                             "sYF", "sSA", "sAD_2to9", "sAD_10up", 
                             "pYF", "pSA", "pAD_2to9", "pAD_10up")) %>% 
      mutate(type = case_when(
        Variable == "BR_all"   ~ "Birth rate",
        Variable == "sPY_all"  ~ "Survival of pouch young",
        Variable == "sYF"      ~ "Survival of young-at-foot",
        Variable == "sSA"      ~ "Survival of subadults",
        Variable == "sAD_2to9" ~ "Survival of adults (2–9)",
        Variable == "sAD_10up" ~ "Survival of adults (10+)",
        Variable == "pYF"      ~ "Proportion of young-at-foot",
        Variable == "pSA"      ~ "Proportion of subadults",
        Variable == "pAD_2to9" ~ "Proportion of adults (2–9)",
        Variable == "pAD_10up" ~ "Proportion of adults (10+)"
      ))
    
    types <- c("Birth rate",
               "Survival of pouch young",
               "Survival of young-at-foot", 
               "Survival of subadults",
               "Survival of adults (2–9)",
               "Survival of adults (10+)",
               "Proportion of young-at-foot",
               "Proportion of subadults", 
               "Proportion of adults (2–9)",
               "Proportion of adults (10+)")
    
    names <- c("Birth\nrate",
               "Survival of\npouch young",
               "Survival of\nyoung-at-foot", 
               "Survival of\nsubadults",
               "Survival of\nadults (2–9)",
               "Survival of\nadults (10+)",
               "Prop. of\nyoung-at-foot",
               "Prop. of\nsubadults", 
               "Prop. of\nadults (2–9)",
               "Prop. of\nadults (10+)")
  }
  
  # order factor levels
  plotData$type <- factor(plotData$type, levels = types)
  plot.colours <- paletteer::paletteer_c("grDevices::Temps", length(unique(plotData$type)))
  
  
  ## Plot contributions --------------------------------------------------------
  
  c.sum <- plotData %>% 
    ggplot(aes(x = type, y = Contribution, group = type)) +
    geom_violin(aes(fill = type, colour = type), alpha = 0.5, scale = "width", draw_quantiles = 0.5) +
    geom_hline(yintercept = 0, colour = "grey70", linetype = "dashed") +
    ylab("Contribution") +
    xlab("") +
    scale_fill_manual(values = plot.colours) +
    scale_colour_manual(values = plot.colours) +
    scale_x_discrete(labels = names) +
    scale_y_continuous(limits = c(-0.008, 0.008),
                       expand = expansion(mult = c(0, 0.02)),
                       breaks = c(-0.008, -0.006, -0.004, -0.002, 0.00, 0.002, 0.004, 0.006, 0.008)) +
    theme_bw() +
    theme(legend.position = "none",
          panel.grid = element_blank(),
          axis.text.x = element_text(size = 10),
          axis.title = element_text(size = 10),
          plot.margin = margin(1, 3, 1, 3))
  
  # ggsave("figures/resultsDave2Covs/LTREsum.jpeg", width = 20.0, height = 12.0, units = c("cm"), dpi = 600)
  
  # survival panel
  S.colours <- c(plot.colours[3:4], rep(plot.colours[5], 8), rep(plot.colours[6], max(0, nAge - 9)))
  names(S.colours) <- c("sYF", "sSA", paste0("sAD_", 2:nAge))
  
  p.S <- ggplot(subset(contData, Variable %in% c("sYF", "sSA", paste0("sAD_", 2:nAge)))) +
    geom_violin(aes(x = factor(Variable, levels = c("sYF", "sSA", paste0("sAD_", 2:nAge))),
                    y = Contribution, fill = Variable), alpha = 0.5, scale = "width", draw_quantiles = 0.5) +
    geom_hline(yintercept = 0, colour = "grey70", linetype = "dashed") +
    ylab("Contribution") +
    xlab("") +
    labs(title = "a) Survival") +
    scale_x_discrete(labels = parse(text = c("italic(S)[0]", "italic(S)[1]", paste0("italic(S)[", 2:nAge, "]")))) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
    scale_fill_manual(values = S.colours) +
    theme_bw() +
    theme(legend.position = "none",
          panel.grid = element_blank(),
          axis.text.x = element_text(size = 10),
          axis.title = element_text(size = 10),
          plot.margin = margin(1, 3, 1, 3))
  
  # birth rate panel
  B.colours <- c(rep(plot.colours[1], nAge - 1))
  names(B.colours) <- c(paste0("BR_", 2:nAge))
  
  p.B <- ggplot(subset(contData, Variable %in% c(paste0("BR_", 2:nAge)))) +
    geom_violin(aes(x = factor(Variable, levels = c(paste0("BR_", 2:nAge))),
                    y = Contribution, fill = Variable), alpha = 0.5, scale = "width", draw_quantiles = 0.5) +
    geom_hline(yintercept = 0, colour = "grey70", linetype = "dashed") +
    ylab("Contribution") +
    xlab("") +
    labs(title = "b) Birth rate") +
    scale_x_discrete(labels = parse(text = paste0("italic(B)[", 2:nAge, "]"))) +
    scale_y_continuous(labels = scales::label_number(),
                       expand = expansion(mult = c(0, 0.02))) +
    scale_fill_manual(values = B.colours) +
    theme_bw() +
    theme(legend.position = "none",
          panel.grid = element_blank(),
          axis.text.x = element_text(size = 10),
          axis.title = element_text(size = 10),
          plot.margin = margin(1, 3, 1, 3))
  
  # reproductive success panel
  R.colours <- c(rep(plot.colours[2], nAge - 1))
  names(R.colours) <- c(paste0("sPY_", 2:nAge))
  
  p.R <- ggplot(subset(contData, Variable %in% c(paste0("sPY_", 2:nAge)))) +
    geom_violin(aes(x = factor(Variable, levels = c(paste0("sPY_", 2:nAge))),
                    y = Contribution, fill = Variable), alpha = 0.5, scale = "width", draw_quantiles = 0.5) +
    geom_hline(yintercept = 0, colour = "grey70", linetype = "dashed") +
    ylab("Contribution") +
    xlab("") +
    labs(title = "c) Survival of pouch young") +
    scale_x_discrete(labels = parse(text = paste0("italic(SP)[", 2:nAge, "]"))) +
    scale_y_continuous(limits = c(0, 0.001), expand = expansion(mult = c(0, 0.02))) +
    scale_fill_manual(values = R.colours) +
    theme_bw() +
    theme(legend.position = "none",
          panel.grid = element_blank(),
          axis.text.x = element_text(size = 10),
          axis.title = element_text(size = 10),
          plot.margin = margin(1, 3, 1, 3))
  
  # population structure panel
  P.colours <- c(rep(plot.colours[length(plot.colours)], nAge + 1))
  names(P.colours) <- c("pYF", "pSA", paste0("pAD_", 2:nAge))
  
  p.P <- ggplot(subset(contData, Variable %in% c("pYF", "pSA", paste0("pAD_", 2:nAge)))) +
    geom_violin(aes(x = factor(Variable, levels = c("pYF", "pSA", paste0("pAD_", 2:nAge))),
                    y = Contribution, fill = Variable), alpha = 0.5, scale = "width", draw_quantiles = 0.5) +
    geom_hline(yintercept = 0, colour = "grey70", linetype = "dashed") +
    ylab("Contribution") +
    xlab("") +
    labs(title = "d) Population proportions") +
    scale_x_discrete(labels = parse(text = c("italic(P)[0]", "italic(P)[1]", paste0("italic(P)[", 2:nAge, "]")))) +
    scale_y_continuous(limits = c(-0.01, 0.01), expand = expansion(mult = c(0, 0.02))) +
    scale_fill_manual(values = P.colours) +
    theme_bw() +
    theme(legend.position = "none",
          panel.grid = element_blank(),
          axis.text.x = element_text(size = 10),
          axis.title = element_text(size = 10),
          plot.margin = margin(1, 3, 1, 3))
  
  
  ## Save plots --------------------------------------------------------------
  
  # e.sum / c.sum
  # ggsave("figures/resultsDave2Covs/elas&ltre.jpeg", width = 18.0, height = 18.0, units = c("cm"), dpi = 600)
  
  pdf(paste0(plotFolder, "/LTRE_sum.pdf"), width = 8, height = 4)
  print(c.sum)
  dev.off()
  
  pdf(paste0(plotFolder, "/LTRE_age.pdf"), width = 8, height = 10)
  print(p.S / p.B / p.R / p.P)
  dev.off()
  
  
  ## Calculate summary -------------------------------------------------------
  
  LTREsummary <- plotData %>%
    group_by(Variable) %>%
    summarise(Mean  = mean(Contribution, na.rm = TRUE),
              Lower = quantile(Contribution, 0.025, na.rm = TRUE),
              Upper = quantile(Contribution, 0.975, na.rm = TRUE),
              .groups = "drop")
  
  
  ## Return --------------------------------------------------------------------
  
  plotList <- c(paste0(plotFolder, "/LTRE_sum.pdf"),
                paste0(plotFolder, "/LTRE_age.pdf"))
  
  if(returnSummary){
    return(c(list(Plots = plotList), list(Summary = LTREsummary)
    ))
  }else{
    return(plotList)
  }
}

