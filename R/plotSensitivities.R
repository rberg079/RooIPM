#' Plot transient sensitivities & elasticities
#'
#' @param sensitivities list. Contains lists of posterior samples of transient sensitivities & elasticities for all parameters.
#' @param nAge integer. Maximum age to consider in the analysis. nAge = 19 by default.
#' @param plotFolder character string. Path to the folder in which to store plots.
#' @param oneProp logical. If TRUE, sums age-structure proportions into a single age structure category. oneProp = TRUE by default.
#' @param returnSummary logical. If TRUE, returns saved file paths & summary data frame. returnSummary = TRUE by default.
#'
#' @returns a character vector of plot names. The plots themselves are saved as pdfs in plotFolder.
#' @export
#'
#' @examples

plotSensitivities <- function(sensitivities, nAge = 19, plotFolder,
                              oneProp = TRUE, returnSummary = TRUE){

  # # for testing purposes
  # sensitivities <- readRDS('results/sensitivities.rds')
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
  
  # make plotting directory if it does not exist
  if(!dir.exists(plotFolder)){
    dir.create(plotFolder)
  }
  
  # function to pull age-specific columns & calculate sums
  processMatrix <- function(mat, prefix) {
    df <- as.data.frame(mat)
    colnames(df) <- paste0(prefix, "_", 1:ncol(df))
    
    df[[paste0(prefix, "_all")]]  <- rowSums(mat)
    df[[paste0(prefix, "_2to9")]] <- rowSums(mat[, 2:9, drop = FALSE])
    df[[paste0(prefix, "_10up")]] <- rowSums(mat[, 10:ncol(mat), drop = FALSE])
    
    return(df)
  }
  
  # initialize empty list for summaries
  summaries <- list()
  
  for(i in 1:2) {
    
    # define names for dynamic labeling & saving
    metric <- ifelse(i == 1, "Sensitivity", "Elasticity")
    
    # select relevant data
    params <- sensitivities[[i]]$samples
    
    # drop sens/elas prefix for generalising
    names(params) <- sub("^.{4}\\.", "", names(params))
    
    
    ## Assemble data -----------------------------------------------------------
    
    # 1D variables
    df_1D <- data.frame(
      sYF = params$sYF,
      sSA = params$sSA,
      pYF = params$pYF,
      pSA = params$pSA
    )
    
    # 2D variables
    df_BR  <- processMatrix(params$BR, "BR")
    df_sPY <- processMatrix(params$sPY, "sPY")
    df_sAD <- processMatrix(params$sAD, "sAD")
    df_pAD <- processMatrix(params$pAD, "pAD")
    
    # bind together
    wideData <- bind_cols(df_1D, df_BR, df_sPY, df_sAD, df_pAD)
    wideData$draw <- 1:nrow(wideData)
    
    # pivot to long format
    sensData <- wideData %>% 
      pivot_longer(-draw, names_to = "Variable", values_to = "Value") %>% 
      mutate(Parameter = str_split_fixed(Variable, "_", 2)[, 1],
             Age = str_split_fixed(Variable, "_", 2)[, 2]) %>% 
      mutate(Age = ifelse(Age == "", NA, Age))
    
    # format for plotting
    if(oneProp){
      tmp <- sensData %>% 
        filter(Variable %in% c("pYF", "pSA", "pAD_2to9", "pAD_10up")) %>% 
        group_by(draw) %>% 
        summarise(Value = sum(Value), .groups = "drop") %>% 
        mutate(Variable = "p_all", type = "Population structure")
      
      plotData <- sensData %>% 
        filter(Variable %in% c("BR_all", "sPY_all", "sYF", "sSA", "sAD_2to9", "sAD_10up")) %>% 
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
                 "Population structure")
      
      names <- c("Birth\nrate",
                 "Survival of\npouch young",
                 "Survival of\nyoung-at-foot", 
                 "Survival of\nsubadults",
                 "Survival of\nadults (2–9)",
                 "Survival of\nadults (10+)",
                 "Population\nstructure")
    }else{
      plotData <- sensData %>% 
        filter(Variable %in% c("BR_all", "sPY_all", "sYF", "sSA", "sAD_2to9", "sAD_10up", 
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
  
  
    ## Plot sensitivities ------------------------------------------------------
    
    plotData$type <- factor(plotData$type, levels = types)
    plot.colours <- paletteer::paletteer_c("grDevices::Temps", length(unique(plotData$type)))
    
    p.sum <- ggplot(plotData, aes(x = type, y = Value, group = type)) +
      geom_violin(aes(fill = type, colour = type), alpha = 0.5, scale = "width", draw_quantiles = 0.5) +
      geom_hline(yintercept = 0, colour = "grey70", linetype = "dashed") +
      ylab(metric) +
      xlab("") +
      scale_fill_manual(values = plot.colours) +
      scale_colour_manual(values = plot.colours) +
      scale_x_discrete(labels = names) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
      theme_bw() +
      theme(legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(size = 10),
            axis.title = element_text(size = 10))
    
    # ggsave("figures/resultsDave2Covs/SENSsum.jpeg", width = 20.0, height = 12.0, units = c("cm"), dpi = 600)
    
    # survival panel
    S.colours <- c(plot.colours[3:4], rep(plot.colours[5], 8), rep(plot.colours[6], max(0, nAge - 9)))
    names(S.colours) <- c("sYF", "sSA", paste0("sAD_", 2:nAge))
    
    p.S <- ggplot(subset(sensData, Variable %in% c("sYF", "sSA", paste0("sAD_", 2:nAge)))) +
      geom_violin(aes(x = factor(Variable, levels = c("sYF", "sSA", paste0("sAD_", 2:nAge))),
                      y = Value, fill = Variable), alpha = 0.5, scale = "width", draw_quantiles = 0.5) +
      geom_hline(yintercept = 0, colour = "grey70", linetype = "dashed") +
      ylab(metric) +
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
    
    p.B <- ggplot(subset(sensData, Variable %in% c(paste0("BR_", 2:nAge)))) +
      geom_violin(aes(x = factor(Variable, levels = c(paste0("BR_", 2:nAge))),
                      y = Value, fill = Variable), alpha = 0.5, scale = "width", draw_quantiles = 0.5) +
      geom_hline(yintercept = 0, colour = "grey70", linetype = "dashed") +
      ylab(metric) +
      xlab("") +
      labs(title = "b) Birth rate") +
      scale_x_discrete(labels = parse(text = paste0("italic(B)[", 2:nAge, "]"))) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
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
    
    p.R <- ggplot(subset(sensData, Variable %in% c(paste0("sPY_", 2:nAge)))) +
      geom_violin(aes(x = factor(Variable, levels = c(paste0("sPY_", 2:nAge))),
                      y = Value, fill = Variable), alpha = 0.5, scale = "width", draw_quantiles = 0.5) +
      geom_hline(yintercept = 0, colour = "grey70", linetype = "dashed") +
      ylab(metric) +
      xlab("") +
      labs(title = "c) Survival of pouch young") +
      scale_x_discrete(labels = parse(text = paste0("italic(SP)[", 2:nAge, "]"))) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
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
    
    p.P <- ggplot(subset(sensData, Variable %in% c("pYF", "pSA", paste0("pAD_", 2:nAge)))) +
      geom_violin(aes(x = factor(Variable, levels = c("pYF", "pSA", paste0("pAD_", 2:nAge))),
                      y = Value, fill = Variable), alpha = 0.5, scale = "width", draw_quantiles = 0.5) +
      geom_hline(yintercept = 0, colour = "grey70", linetype = "dashed") +
      ylab(metric) +
      xlab("") +
      labs(title = "d) Population proportions") +
      scale_x_discrete(labels = parse(text = c("italic(P)[0]", "italic(P)[1]", paste0("italic(P)[", 2:nAge, "]")))) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
      scale_fill_manual(values = P.colours) +
      theme_bw() +
      theme(legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(size = 10),
            axis.title = element_text(size = 10),
            plot.margin = margin(1, 3, 1, 3))
    
    
    ## Save plots --------------------------------------------------------------
    
    # p.S / p.B / p.R / p.P
    # ggsave("figures/resultsDave2Covs/SENSage2.jpeg", width = 20.0, height = 24.0, units = c("cm"), dpi = 600)
    
    pdf(paste0(plotFolder, ifelse(i == 1, "/Sensitivities", "/Elasticities"), "_sum.pdf"), width = 8, height = 4)
    print(p.sum)
    dev.off()
    
    pdf(paste0(plotFolder, ifelse(i == 1, "/Sensitivities", "/Elasticities"), "_age.pdf"), width = 8, height = 10)
    print(p.S / p.B / p.R / p.P)
    dev.off()
    
    
    ## Calculate summary -------------------------------------------------------
    
    sensSummary <- sensData %>%
      group_by(Variable) %>%
      summarise(Mean  = mean(Value, na.rm = TRUE),
                Lower = quantile(Value, 0.025, na.rm = TRUE),
                Upper = quantile(Value, 0.975, na.rm = TRUE),
                .groups = "drop")
    
    summaries[[metric]] <- sensSummary
  }
  
  
  ## Return --------------------------------------------------------------------
  
  plotList <- c(paste0(plotFolder, c("/Sensitivities", "/Elasticities"), "_sum.pdf"),
                paste0(plotFolder, c("/Sensitivities", "/Elasticities"), "_age.pdf"))
  
  # Return plots alone, or seamlessly tack the summaries on to the same level
  if(returnSummary){
    return(c(list(Plots = plotList), summaries))
  }else{
    return(plotList)
  }
}
  
  