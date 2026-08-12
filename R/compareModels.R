#' Compare outputs of different models
#'
#' @param nYear integer. Number of years to consider in the analysis. nYear = 18 by default.
#' @param minYear integer. First year to consider in the analysis. minYear = 2008 by default.
#' @param maxYear integer. Last year to consider in the analysis. maxYear = minYear + nYear - 1 by default.
#' @param nAgeC.S integer. Number of age classes to consider in the survival model. nAgeC.S = 12 by default.
#' @param plotAges integer vector. Ages to plot in time series plots. plotAges = c(2, 6, 10, 14) by default.
#' @param plotYears integer vector. Years to plot in density plots. plotYears = c(2, 6, 10, 14) by default.
#' @param postPaths character vector. Paths to .rds files containing posterior samples from models to compare.
#' @param modelNames character vector. User-defined names for models to compare. 
#' @param plotFolder character string. Path to the folder in which to store plots.
#' @param returnSumData logical. If TRUE, returns a data frame containing posterior samples from all 
#' compared models as an object in the R global environment. If FALSE (default), no data is returned.
#'
#' @returns a data frame of posterior samples from all compared models, provided returnSumData is TRUE.
#' @export
#'
#' @examples

compareModels <- function(nYear = 18, minYear = 2008, maxYear, nAgeC.S = 12,
                          plotAges = c(2, 6, 10, 14), plotYears = c(2, 6, 10, 14),
                          postPaths, modelNames, plotFolder, returnSumData = FALSE){
  
  # # for testing purposes
  # nYear = 18
  # minYear = 2008
  # maxYear = minYear + nYear - 1
  # nAgeC.S = 12
  # plotAges = c(2, 6, 10, 14)
  # plotYears = c(2, 6, 10, 14)
  # postPaths = c("results/IPM_CJSen_RSen_AB_DynDens_dCJS_12_noW_25BR_Dave2Covs_stochV_8chains.rds",
  #               "results/IPM_CJSen_RSen_AB_DynDens_dCJS_12_noW_25BR_Dave3Covs3REs_stochV_8chains.rds"
  #               )
  # modelNames = c("IPM_Dave2covs",
  #                "IPM_Dave3covs3REs"
  #                )
  # plotFolder = c("figures/densityChecks/final2")
  # returnSumData = TRUE
  # nModels <- length(modelNames)
  
  
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
  
  # set maxYear if not provided
  if(missing(maxYear)){
    maxYear <- minYear + nYear - 1
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
                               grepl('EpsilonT|nYF|nSA|nTOT|sYF|sSA|propF|dens.true|veg.true', Parameter) ~ Idx1,
                               grepl('nAD|BR|sPY|sAD|S', Parameter) ~ Idx2),
           AgeIdx  = case_when(grepl('Beta|EpsilonI|EpsilonT|SigmaT|Mu.O|nYF|nSA|nTOT|sYF|sSA|propF', Parameter) ~ NA_real_,
                               grepl('Mu.S|Mu.B|Mu.R|nAD|BR|sPY|sAD|S', Parameter) ~ Idx1),
           Year = YearIdx + minYear - 1,
           Age  = case_when(grepl('Mu.S|Mu.B|Mu.R|nAD|BR|sPY|sAD|S', Parameter) ~ AgeIdx,
                            grepl('nYF|sYF', Parameter) ~ 0,
                            grepl('nSA|sSA', Parameter) ~ 1,
                            TRUE ~ NA_real_),
           ParamName = word(Parameter, 1, sep = "\\["),
           ParamName = ifelse(ParamName %in% c('nAD', 'sAD', 'BR', 'sPY') & AgeIdx %in% plotAges,
                              paste0(ParamName, '[', AgeIdx, ']'), ParamName))
  
  sum.dat <- sum.dat %>%
    left_join(idx.dat, by = "Parameter")
  
  
  ## Prep parameters for plotting ----------------------------------------------
  
  # set parameter groups for plotting posterior density overlaps
  plot.params <- list(
    CJS_mus = c(paste0('Mu.S[', 1:nAgeC.S, ']')),
    
    CJS_covs = c('BetaD.S', 'BetaV.S',
                 'BetaD.Sy', 'BetaD.Sp', 'BetaD.So',
                 'BetaV.Sy', 'BetaV.Sp', 'BetaV.So'),
    
    CJS_REs = c(paste0('EpsilonT.S[', plotYears, ']'),
                'SigmaT.S', 'SigmaT.Sy', 'SigmaT.Sp', 'SigmaT.So'),
    
    CJS_obs = c('Mu.O', 'SigmaT.O', paste0('EpsilonT.O[', 1:nYear, ']')),
    
    RS_muB = c(paste0('Mu.B[', 1:nAgeC.S, ']')),
    
    RS_BR = c(expand.grid(a = plotAges, t = plotYears) %>%
                mutate(param = paste0('BR[', a, ', ', t, ']')) %>%
                pull(param)),
    
    RS_muR = c(paste0('Mu.R[', 1:nAgeC.S, ']')),
    
    RS_Ra = c(expand.grid(a = plotAges, t = plotYears) %>%
                mutate(param = paste0('sPY[', a, ', ', t, ']')) %>%
                pull(param)),
    
    RS_covs = c('BetaD.R', 'BetaD.Rp', 'BetaD.Ro'),
    
    RS_REs = c(paste0('EpsilonT.R[', plotYears, ']'),
               paste0('EpsilonT.B[', plotYears, ']'),
               'SigmaI.R', 'SigmaT.B', 'SigmaT.R',
               'SigmaT.Rp', 'SigmaT.Ro'),
    
    POP_NAs = c(expand.grid(a = plotAges, t = plotYears) %>%
                  mutate(param = paste0('nAD[', a, ', ', t, ']')) %>%
                  pull(param)),
    
    POP_NTs = c(paste0('nYF[', plotYears, ']'),
                paste0('nSA[', plotYears, ']'),
                paste0('nTOT[', plotYears, ']')))
  
  # set parameters for plotting time series of posterior summaries
  plotTS.VRs <- list(
    ParamNames = c(expand.grid(a = plotAges) %>% 
                     mutate(param = paste0('BR[', a, ']')) %>%
                     pull(param),
                   expand.grid(a = plotAges) %>% 
                     mutate(param = paste0('sPY[', a, ']')) %>%
                     pull(param),
                   'sYF',
                   'sSA', 
                   expand.grid(a = plotAges) %>% 
                     mutate(param = paste0('sAD[', a, ']')) %>%
                     pull(param),
                   'dens.true',
                   'veg.true'),
    
    ParamLabels = c(expand.grid(a = plotAges) %>% 
                      mutate(name = paste0('Breeding rate (age ', a, ')')) %>% 
                      pull(name),
                    expand.grid(a = plotAges) %>% 
                      mutate(name = paste0('Survival to pouch exit (', a, ' y/o moms)')) %>% 
                      pull(name),
                    'Survival of young-at-foot (0 yrs)',
                    'Survival of subadults (1 yr)',
                    expand.grid(a = plotAges) %>% 
                      mutate(name = paste0('Survival of adults(', a, ' yrs)')) %>% 
                      pull(name),
                    '"True" density',
                    '"True" forage'))
  
  plotTS.Ns <- list(
    ParamNames = c('nYF',
                   'nSA',
                   expand.grid(a = plotAges) %>% 
                     mutate(param = paste0('nAD[', a, ']')) %>% 
                     pull(param),
                   'nTOT'),
    
    ParamLabels = c('# of young-at-foot',
                    '# of subadults',
                    expand.grid(a = plotAges) %>% 
                      mutate(name = paste0('# of adults(', a, ' yrs)')) %>% 
                      pull(name),
                    'Total # of females'))
  
  # set plotting colors
  plot.cols <- paletteer_c("grDevices::Temps", nModels)
  
  
  ## Plot ----------------------------------------------------------------------
  
  # posterior overlaps
  pdf(paste0(plotFolder, "/PostDensities.pdf"), width = 9, height = 6)
  for(x in 1:length(plot.params)){
    
    print(
      ggplot(subset(post.dat, Parameter %in% plot.params[[x]]),
             aes(x = Value, color = Model, fill = Model)) + 
        geom_density(alpha = 1/nModels) + 
        facet_wrap(~Parameter, scales = "free") + 
        scale_fill_manual(values = plot.cols) +
        scale_color_manual(values = plot.cols) + 
        theme_bw() + theme(panel.grid = element_blank())
    )
  }
  dev.off()
  
  # posterior summary time series of vital rates
  pdf(paste0(plotFolder, "/PostSummariesTS_VRs.pdf"), width = 8, height = 4)
  for(x in 1:length(plotTS.VRs$ParamNames)){
    
    print(
      ggplot(subset(sum.dat, ParamName == plotTS.VRs$ParamNames[x] & Year >= minYear & Year <= maxYear), aes(group = Model)) + 
        geom_ribbon(aes(x = Year, ymin = Lower, ymax = Upper, fill = Model), alpha = 1/nModels) + 
        geom_line(aes(x = Year, y = Median, color = Model)) + 
        scale_fill_manual(values = plot.cols) +
        scale_color_manual(values = plot.cols) + 
        scale_x_continuous(breaks = c(minYear:maxYear), labels = c(minYear:maxYear)) + 
        ggtitle(plotTS.VRs$ParamLabels[x]) +  
        theme_bw() + theme(panel.grid.minor = element_blank(), 
                           panel.grid.major.y = element_blank(), 
                           axis.text.x = element_text(angle = 45, vjust = 0.5))
    )
    
  }
  dev.off()
  
  # posterior summary time series of population sizes
  pdf(paste0(plotFolder, "/PostSummariesTS_Ns.pdf"), width = 8, height = 4)
  for(x in 1:length(plotTS.Ns$ParamNames)){
    
    print(
      ggplot(subset(sum.dat, ParamName == plotTS.Ns$ParamNames[x] & Year > minYear & Year <= maxYear), aes(group = Model)) + 
        geom_ribbon(aes(x = Year, ymin = Lower, ymax = Upper, fill = Model), alpha = 1/nModels) + 
        geom_line(aes(x = Year, y = Median, color = Model)) + 
        scale_fill_manual(values = plot.cols) +
        scale_color_manual(values = plot.cols) + 
        scale_x_continuous(breaks = c(minYear:maxYear), labels = c(minYear:maxYear)) + 
        ggtitle(plotTS.Ns$ParamLabels[x]) +  
        theme_bw() + theme(panel.grid.minor = element_blank(), 
                           panel.grid.major.y = element_blank(), 
                           axis.text.x = element_text(angle = 45, vjust = 0.5))
    )
    
  }
  dev.off()
  
  # optional: return summary data
  if(returnSumData){
    return(sum.dat)
  }
}

