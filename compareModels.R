#' Compare outputs of different models
#'
#' @param nAge integer. Maximum age to consider in the analysis. nAge = 22 by default.
#' @param nAgeC integer. Number of age classes to consider in the analysis. nAgeC = 5 by default.
#' @param nYear integer. Number of years to consider in the analysis. nYear = 17 by default.
#' @param minYear integer. First year to consider in the analysis. minYear = 2008 by default.
#' @param maxYear integer. Last year to consider in the analysis. maxYear = minYear + nYear - 1 by default.
#' @param postPaths character vector. Paths to .rds files containing posterior samples from models to compare.
#' @param modelNames character vector. User-defined names for models to compare. 
#' @param plotFolder character string. Path to the folder in which to store comparison plots.
#' @param returnSumData logical. If TRUE, returns a data frame containing posterior samples from all 
#' compared models as an object in the R global environment. If FALSE (default), no data is returned.
#'
#' @returns a data frame of posterior samples from all compared models, provided returnSumData is TRUE.
#' @export
#'
#' @examples

compareModels <- function(nAge = 19, nAgeC = 5, nYear = 17, nNoProp = 5, minYear = 2008, maxYear,
                          postPaths, modelNames, plotFolder, returnSumData = FALSE){
  
  # # for testing purposes
  # nAge = 19
  # nAgeC = 5
  # nYear = 17
  # nNoProp = 5
  # minYear = 2008
  # maxYear = minYear + nYear - 1
  # postPaths = c("results/IPM_CJSen_RSen.rds", "results/IPM_CJSen_RSen_AB.rds")
  # modelNames = c("CJSen/RSen", "CJSen/RSen/AB")
  # plotFolder = c("figures/AB")
  # returnSumData = TRUE
  # nModels <- length(modelNames)

  ## Set up --------------------------------------------------------------------

  library(coda)
  suppressPackageStartupMessages(library(tidyverse))
  suppressPackageStartupMessages(library(data.table))
  library(NatParksPalettes)
  
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
           YearIdx = case_when(grepl('^Sigma\\.S\\[|^Xi\\.S\\[|^EpsilonI\\.Ri\\[|^Mu\\.Ri\\[|^Mu\\.Ra\\[', Parameter) ~ NA_real_,
                               !is.na(Idx2) & !grepl('^Gamma\\.S\\[', Parameter) ~ Idx2,
                               !is.na(Idx1) & grepl('^Gamma\\.S\\[', Parameter) ~ Idx1,
                               !is.na(Idx1) & !grepl('Beta', Parameter) ~ Idx1),
           AgeIdx = case_when(grepl('^Sigma\\.S\\[|^Xi\\.S\\[', Parameter) ~ Idx1,
                              !is.na(Idx2) & !grepl('^Gamma\\.S\\[', Parameter) ~ Idx1,
                              !is.na(Idx2) & grepl('^Gamma\\.S\\[', Parameter) ~ Idx2,
                              is.na(Idx2) & grepl('Beta|Mu.Ri|Mu.Ra', Parameter) ~ Idx1),
           Year = YearIdx + minYear - 1,
           Age = case_when(grepl('^nSA\\[|^nAD\\[|^Ra\\[|^sSA\\[|^sAD\\[', Parameter) ~ AgeIdx,
                           grepl('^nSA\\[|^sSA\\[', Parameter) ~ 1,
                           grepl('^nYAF\\[|^sYAF\\[', Parameter) ~ 0,
                           TRUE ~ NA_real_),
           # AgeClass = ifelse(grepl('^BetaA\\.S\\[|^BetaD\\.S\\[|^BetaV\\.S\\[|^S\\[', Parameter), AgeIdx, NA),
           AgeClass = ifelse(grepl('Beta|^S\\[', Parameter), AgeIdx, NA),
           ParamName = word(Parameter, 1, sep = "\\["))
  
  sum.dat <- sum.dat %>%
    left_join(idx.dat, by = "Parameter")
  
  
  ## Prep parameters for plotting ----------------------------------------------
  
  # set parameter groups for plotting posterior density overlaps
  plot.params <- list(
    CJScovEF = c(paste0('BetaA.S[', 1:nAgeC, ']'),
                 paste0('BetaD.S[', 1:nAgeC, ']'),
                 paste0('BetaV.S[', 1:nAgeC, ']')),
    
    CJSranEF = c(paste0('Sigma.S[', 1:nAgeC, ', ', 1:nAgeC, ']')),
    
    CJSestO = c('Mu.O', 'Sigma.O', paste0('Epsilon.O[', 1:nYear, ']')),
    
    CJSestS = c(expand.grid(a = 1:nAgeC, t = c(2, 6, 10, 14)) %>% 
                  mutate(param = paste0('S[', a, ', ', t, ']')) %>% 
                  pull(param)),
    
    RScovEF = c('BetaD.R', 'BetaV.R', 'BetaW.R'),

    RSranEF = c(paste0('EpsilonT.Ri[', c(2, 6, 10, 14), ']'),
                paste0('EpsilonT.Ra[', c(2, 6, 10, 14), ']'),
                paste0('EpsilonT.B[', c(2, 6, 10, 14), ']'),
                'SigmaT.Ri', 'SigmaT.Ra', 'SigmaT.B'),
    
    RSestBt = c(paste0('Bt[', 1:(nYear-1), ']')),
    
    RSestRa = c(expand.grid(a = c(2, 6, 10, 14), t = c(2, 6, 10, 14)) %>%
                  mutate(param = paste0('Ra[', a, ', ', t, ']')) %>%
                  pull(param)),
    
    POPestNA = c(expand.grid(a = c(2, 6, 10, 14), t = c(2, 6, 10, 14)) %>%
                   mutate(param = paste0('nAD[', a, ', ', t, ']')) %>%
                   pull(param)),
    
    POPestNT = c(paste0('nYAF[', c(2, 6, 10, 14), ']'),
                 paste0('nSA[', c(2, 6, 10, 14), ']'),
                 paste0('nTOT[', c(2, 6, 10, 14), ']')),
    
    ABestAB = c(paste0('ab[', 1:nYear, ']')))
  
  # set parameters for plotting time series of posterior summaries
  plotTS.VRs <- list(
    ParamNames = c('Bt', 'Ra', 'sYAF', 'sSA', 'sAD'),
    ParamLabels = c('Breeding rate',
                    'Survival to pouch exit',
                    'Survival of young-at-foot (0 yrs)',
                    'Survival of subadults (1 yr)',
                    'Survival of adults (2+ yrs)'))
  
  plotTS.Ns <- list(
    ParamNames = c('nYAF', 'nSA', 'nAD', 'nTOT', 'ab'),
    ParamLabels = c('# of young-at-foot', '# of subadults',
                    '# of adults', 'Total # of females', 'Abundance'))
  
  # set plotting colors
  # plot.cols <- paletteer_c("grDevices::Temps", length(modelNames))
  plot.cols <- natparks.pals("Banff", nModels)
  
  
  ## Plot ----------------------------------------------------------------------
  
  # posterior overlaps
  pdf(paste0(plotFolder, "/PostDensities.pdf"), width = 9, height = 6)
  for(x in 1:length(plot.params)){
    
    print(
      ggplot(subset(post.dat, Parameter %in% plot.params[[x]]),
             aes(x = Value, color = Model, fill = Model)) + 
        geom_density(alpha = 1/nModels) + 
        facet_wrap(~Parameter, scales = "free") + 
        # scale_fill_viridis_d() +
        # scale_color_viridis_d() + 
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

