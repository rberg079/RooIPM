#' Plot age-specific vital rate intercepts
#'
#' @param MCMCsamples mcmc.list or matrix containing the output of the fitted IPM.
#' @param nAge integer. Number of ages, or maximum age, to consider. nAge = 19 by default.
#' @param plotFolder character string. Path to the folder in which to store plots.
#' @param returnSummary logical. If TRUE, returns saved file paths and summary data frame. returnSummary = TRUE by default.
#'
#' @return a character vector of plot names & optionally the intercept summary metrics. The plots themselves are saved as pdfs in plotFolder.
#' @export
#'
#' @examples

plotIPM_VitalRateIntercepts <- function(MCMCsamples, nAge = 19, plotFolder, returnSummary = TRUE){
  
  # # for testing purposes
  # MCMCsamples <- readRDS('results/IPM_CJSen_RSen_AB_DynDens_dCJS_12_noW_25BR_S33_R22_stochV_8chains.rds')
  # plotFolder = "figures/densityChecks/varyNvarsSR"
  # nAge = 19
  # returnSummary <- TRUE
  
  
  ## Set up --------------------------------------------------------------------
  
  suppressPackageStartupMessages(library(tidyverse))
  suppressPackageStartupMessages(library(patchwork))
  library(paletteer)
  library(scales)
  library(coda)
  
  # make plotting directory if it does not exist already
  if(!dir.exists(plotFolder)){
    dir.create(plotFolder, recursive = TRUE)
  }
  
  # convert MCMC samples to matrix
  out.mat <- as.matrix(MCMCsamples)
  
  # color palette matching LTRE style
  plot.colours <- paletteer::paletteer_c("grDevices::Temps", 3)
  col.S <- plot.colours[1]
  col.B <- plot.colours[2]
  col.R <- plot.colours[3]
  
  
  ## Extract age-indexed parameters --------------------------------------------
  
  get_param_df <- function(param_prefix, ages){
    res_list <- list()
    for(a in ages){
      pname <- paste0(param_prefix, "[", a, "]")
      if(pname %in% colnames(out.mat)){
        res_list[[pname]] <- data.frame(
          Draw = 1:nrow(out.mat),
          Age = a,
          Variable = pname,
          Value = out.mat[, pname]
        )
      }
    }
    if(length(res_list) > 0){
      return(bind_rows(res_list))
    }else{
      return(NULL)
    }
  }
  
  # extract data across specified age range
  df.S <- get_param_df("Mu.S", 1:nAge)
  df.B <- get_param_df("Mu.B", 1:nAge)
  df.R <- get_param_df("Mu.R", 1:nAge)
  
  plot_list <- list()
  
  
  ## Panel a) Survival (Mu.S) --------------------------------------------------
  
  if(!is.null(df.S)){
    s_ages <- unique(df.S$Age-1)
    df.S$Variable <- factor(df.S$Variable, levels = paste0("Mu.S[", s_ages, "]"))
    
    p.S <- ggplot(df.S, aes(x = Variable, y = Value)) +
      geom_violin(fill = col.S, colour = col.S, alpha = 0.5, scale = "width", draw_quantiles = 0.5) +
      ylab("Baseline probability") +
      xlab("") +
      labs(title = "a) Survival intercept (Mu.S)") +
      scale_x_discrete(labels = parse(text = paste0("italic(S)[", s_ages, "]"))) +
      scale_y_continuous(labels = scales::label_number(), expand = expansion(mult = c(0.02, 0.02))) +
      theme_bw() +
      theme(panel.grid = element_blank(),
            axis.text.x = element_text(size = 10),
            axis.title = element_text(size = 10),
            plot.margin = margin(1, 3, 1, 3))
    
    plot_list[["p.S"]] <- p.S
  }
  
  
  ## Panel b) Birth rate (Mu.B) ------------------------------------------------
  
  if(!is.null(df.B)){
    b_ages <- unique(df.B$Age+1)
    df.B$Variable <- factor(df.B$Variable, levels = paste0("Mu.B[", b_ages, "]"))
    
    p.B <- ggplot(df.B, aes(x = Variable, y = Value)) +
      geom_violin(fill = col.B, colour = col.B, alpha = 0.5, scale = "width", draw_quantiles = 0.5) +
      ylab("Baseline probability") +
      xlab("") +
      labs(title = "b) Birth rate intercept (Mu.B)") +
      scale_x_discrete(labels = parse(text = paste0("italic(B)[", b_ages, "]"))) +
      scale_y_continuous(labels = scales::label_number(), expand = expansion(mult = c(0.02, 0.02))) +
      theme_bw() +
      theme(panel.grid = element_blank(),
            axis.text.x = element_text(size = 10),
            axis.title = element_text(size = 10),
            plot.margin = margin(1, 3, 1, 3))
    
    plot_list[["p.B"]] <- p.B
  }
  
  
  ## Panel c) Survival of pouch young (Mu.R) -----------------------------------
  
  if(!is.null(df.R)){
    r_ages <- unique(df.R$Age+1)
    df.R$Variable <- factor(df.R$Variable, levels = paste0("Mu.R[", r_ages, "]"))
    
    p.R <- ggplot(df.R, aes(x = Variable, y = Value)) +
      geom_violin(fill = col.R, colour = col.R, alpha = 0.5, scale = "width", draw_quantiles = 0.5) +
      ylab("Baseline probability") +
      xlab("") +
      labs(title = "c) Survival of pouch young intercept (Mu.R)") +
      scale_x_discrete(labels = parse(text = paste0("italic(sPY)[", r_ages, "]"))) +
      scale_y_continuous(labels = scales::label_number(), expand = expansion(mult = c(0.02, 0.02))) +
      theme_bw() +
      theme(panel.grid = element_blank(),
            axis.text.x = element_text(size = 10),
            axis.title = element_text(size = 10),
            plot.margin = margin(1, 3, 1, 3))
    
    plot_list[["p.R"]] <- p.R
  }
  
  
  ## Combine & Save Plot -------------------------------------------------------
  
  combined_plot <- patchwork::wrap_plots(plot_list, ncol = 1)
  
  out_path <- paste0(plotFolder, "/VRintercepts.pdf")
  pdf(out_path, width = 8, height = 4 * length(plot_list))
  print(combined_plot)
  dev.off()
  
  
  ## Calculate Summary ---------------------------------------------------------
  
  all_data <- bind_rows(df.S, df.B, df.R)
  
  interceptSummary <- all_data %>%
    group_by(Variable) %>%
    summarise(Mean   = mean(Value, na.rm = TRUE),
              Median = median(Value, na.rm = TRUE),
              Lower  = quantile(Value, 0.025, na.rm = TRUE),
              Upper  = quantile(Value, 0.975, na.rm = TRUE),
              .groups = "drop")
  
  
  ## Return --------------------------------------------------------------------
  
  if(returnSummary){
    return(list(Plots = out_path, Summary = interceptSummary))
  }else{
    return(out_path)
  }
}

