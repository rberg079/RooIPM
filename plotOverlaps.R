#' Assess prior and posterior estimate overlap
#'
#' @param nYear integer. Number of years to consider in the analysis. nYear = 17 by default.
#' @param minYear integer. First year to consider in the analysis. minYear = 2008 by default.
#' @param nAgeC.S integer. Number of age classes to consider in the survival model. nAgeC.S = 6 by default.
#' @param nAgeC.R integer. Number of age classes to consider in the reproductive success model. nAgeC.R = 6 by default.
#' @param plotAges integer vector. Ages to plot in time series plots. plotAges = c(2, 6, 10, 14) by default.
#' @param plotYears integer vector. Years to plot in density plots. plotYears = c(2, 6, 10, 14) by default.
#' @param postPaths character vector. Paths to .rds files containing posterior samples from models to compare.
#' @param modelNames character vector. User-defined names for models to compare. 
#' @param plotFolder character string. Path to the folder in which to store plots.
#'
#' @returns a pdf containing a series of graphs of prior and posterior overlaps.
#' @export
#'
#' @examples

plotOverlaps <- function(nYear = 17, minYear = 2008, nAgeC.S = 6, nAgeC.R = 6,
                         plotAges = c(2, 6, 10, 14), plotYears = c(2, 6, 10, 14),
                         postPaths, modelNames, plotFolder){
  
  # for testing purposes
  nYear = 17
  minYear = 2008
  nAgeC.S = 6
  nAgeC.R = 6
  plotAges = c(2, 6, 10, 14)
  plotYears = c(2, 6, 10, 14)
  postPaths = c("results/IPM_CJSen_RSen_AB_DynDens_dCJS.rds")
  modelNames = c("IPM_dCJS")
  plotFolder = c("figures")
  nModels = 1
  
  
  ## Set up --------------------------------------------------------------------
  
  library(coda)
  suppressPackageStartupMessages(library(tidyverse))
  suppressPackageStartupMessages(library(data.table))
  library(NatParksPalettes)
  library(paletteer)
  
  # check that models are specified correctly
  if(missing(postPaths)){
    stop("Models must be specified via file paths (postPaths).")
  }
  
  # make plotting directory if it does not exist already
  if(!dir.exists(plotFolder)){
    dir.create(plotFolder)
  }
  
  # count number of models
  nModels <- length(modelNames)
  
  
  ## Reformat posterior samples ------------------------------------------------
  
  post.dat <- data.frame()
  for(i in 1:nModels){
    
    # read in RDS as mcmc.list
    post <- readRDS(postPaths[i])
    
    # convert to matrix then data.table
    samples <- do.call(rbind, lapply(post, as.matrix))
    rownames(samples) <- 1:nrow(samples)  # add row names
    model.dat <- as.data.table(samples, keep.rownames = "Sample")
    
    # reshape to long format
    model.dat <- melt(model.dat,
                      id.vars = "Sample",
                      variable.name = "Parameter",
                      value.name = "Value")
    
    # add identifier & bind
    model.dat[, Model := modelNames[i]]
    post.dat <- rbindlist(list(post.dat, model.dat))
  }
  
  # summarize samples as median & 95% CI
  sum.dat <- post.dat %>%
    group_by(Parameter, Model) %>%
    summarise(Median = median(Value, na.rm = TRUE),
              Lower = quantile(Value, probs = 0.025, na.rm = TRUE),
              Upper = quantile(Value, probs = 0.975, na.rm = TRUE),
              .groups = "keep") %>%
    mutate(Parameter = as.character(Parameter)) %>%
    ungroup()
  
  # extract age & year information
  idx.dat <- data.frame(cbind(unique(sum.dat$Parameter),
                              str_extract_all(unique(sum.dat$Parameter),
                                              pattern = "\\d+", simplify = TRUE))) %>%
    rename("Parameter" = "X1", "Idx1" = "X2", "Idx2" = "X3") %>%
    mutate(Idx1 = as.numeric(ifelse(Idx1 %in% c("", 0), NA, Idx1)),
           Idx2 = as.numeric(ifelse(Idx2 %in% c("", 0), NA, Idx2)),
           YearIdx = case_when(grepl('Beta|EpsilonI|Mu|Sigma', Parameter) ~ NA_real_,
                               grepl('Bt|EpsilonT|nYF|nSA|nTOT|propF|sYF|sSA', Parameter) ~ Idx1,
                               grepl('nAD|sPY|S|sAD', Parameter) ~ Idx2),
           AgeIdx  = case_when(grepl('BetaD.R|BetaV.R|BetaW.R|Bt|EpsilonI|EpsilonT|Mu.B|Mu.O|nYF|nSA|nTOT|propF|SigmaT|sYF|sSA', Parameter) ~ NA_real_,
                               grepl('BetaD.S|BetaV.S|BetaW.S|Mu.S|Mu.R|nAD|S|sPY|sAD', Parameter) ~ Idx1),
           Year = YearIdx + minYear - 1,
           Age  = case_when(grepl('Mu.R|nAD|sPY|sAD', Parameter) ~ AgeIdx,
                            grepl('nYF|sYF', Parameter) ~ 0,
                            grepl('nSA|sSA', Parameter) ~ 1,
                            TRUE ~ NA_real_),
           ParamName = word(Parameter, 1, sep = "\\["),
           ParamName = ifelse(ParamName %in% c('nAD', 'sAD', 'sPY') & AgeIdx %in% plotAges,
                              paste0(ParamName, '[', AgeIdx, ']'), ParamName))
  
  sum.dat <- sum.dat %>%
    left_join(idx.dat, by = "Parameter")
  
  
  ## Generate prior samples ----------------------------------------------------
  
  niter <- nrow(samples) 
  
  # create a list of prior draws
  prior_list <- list()
  
  # survival priors
  for(a in 1:nAgeC.S) prior_list[[paste0("Mu.S[", a, "]")]] <- runif(niter, 0, 1)
  prior_list[["BetaD.S"]] <- runif(niter, -5, 5)
  prior_list[["BetaV.S"]] <- runif(niter, -5, 5)
  prior_list[["BetaW.S"]] <- runif(niter, -5, 5)
  prior_list[["SigmaT.S"]] <- runif(niter, 0, 10)
  
  prior_list[["Mu.O"]] <- runif(niter, 0.01, 0.99)
  prior_list[["SigmaT.O"]] <- runif(niter, 0.01, 10)
  
  # reproductive success priors
  for(a in 1:nAgeC.R) prior_list[[paste0("Mu.R[", a, "]")]] <- runif(niter, 0, 1)
  prior_list[["Mu.B"]] <- runif(niter, 0, 1)
  prior_list[["BetaD.R"]] <- runif(niter, -5, 5)
  prior_list[["SigmaI.R"]] <- runif(niter, 0, 10)
  prior_list[["SigmaT.R"]] <- runif(niter, 0, 10)
  prior_list[["SigmaT.B"]] <- runif(niter, 0, 10)
  
  # convert to data.table
  # reshape to long format (matching post.dat)
  prior.df <- as.data.table(prior_list)
  prior.df[, Sample := 1:.N]
  prior.dat <- melt(prior.df, 
                    id.vars = "Sample", 
                    variable.name = "Parameter", 
                    value.name = "Value")
  
  # tag as the "Prior" model
  prior.dat[, Model := "Prior"]
  
  # bind to the posterior data
  overlap.dat <- rbindlist(list(post.dat, prior.dat), fill = TRUE)
  
  
  ## Prep parameters for plotting ----------------------------------------------
  
  # set parameter groups for plotting overlaps
  plot.params <- list(
    
    CJS = c(paste0('Mu.S[', 1:nAgeC.S, ']'),
            'BetaD.S', 'BetaV.S', 'BetaW.S',
            'SigmaT.S', 'Mu.O', 'SigmaT.O'),
    
    RS  = c('BetaD.R', 'BetaV.R', 'BetaW.R',
            'SigmaT.R', 'SigmaT.B'))
  
  # set plotting colors
  plot.cols <- paletteer_c("grDevices::Temps", nModels)
  
  
  ## Plot ----------------------------------------------------------------------
  
  # prior & posterior overlaps
  overlap <- c("Prior", modelNames)
  plot.cols.overlap <- c("Prior" = "grey80", 
                         setNames(paletteer_c("grDevices::Temps", nModels), modelNames))
  
  pdf(paste0(plotFolder, "/PriorPostOverlaps.pdf"), width = 9, height = 6)
  
  for(x in 1:length(plot.params)){
    
    plot_subset <- subset(overlap.dat, Parameter %in% plot.params[[x]])
    
    if(nrow(plot_subset) > 0) {
      print(
        ggplot(plot_subset, aes(x = Value, color = Model, fill = Model)) +
          geom_density(alpha = 0.4) + 
          facet_wrap(~Parameter, scales = "free") + 
          scale_fill_manual(values = plot.cols.overlap) +
          scale_color_manual(values = plot.cols.overlap) + 
          theme_bw() + 
          theme(panel.grid = element_blank())
      )
    }
  }
  dev.off()
  
}

