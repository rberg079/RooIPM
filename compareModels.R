#' Compare outputs of different models
#'
#' @param nAge integer. Maximum age to consider in the analysis. nAge = 22 by default.
#' @param nAgeC integer. Number of age classes to consider in the analysis. nAgeC = 5 by default.
#' @param ntimes integer. Number of years to consider in the analysis. ntimes = 17 by default.
#' @param minYear integer. First year to consider in the analysis. minYear = 2008 by default.
#' @param maxYear integer. Last year to consider in the analysis. maxYear = minYear + ntimes - 1 by default.
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

compareModels <- function(nAge = 22, nAgeC = 5, ntimes = 17, minYear = 2008, maxYear,
                          postPaths, modelNames, plotFolder, returnSumData = FALSE){
  

  ## Set up --------------------------------------------------------------------

  library(coda)
  library(tidyverse)
  library(data.table)
  # library(paletteer)
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
    maxYear <- minYear + ntimes - 1
  }
  
  # count number of models
  nModels <- length(modelNames)
  
  
  ## Reformat posterior samples ------------------------------------------------
  
  post.dat <- data.frame()
  for(i in 1:nModels){
    
    # read in RDS as mcmc.list
    post <- readRDS(postPaths[i])$out.mcmc
    
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
           YearIdx = case_when(grepl('^Sigma\\.raw\\[|^xi\\[', Parameter) ~ NA_real_,
                               !is.na(Idx2) & !grepl('gamma', Parameter) ~ Idx2,
                               !is.na(Idx1) & grepl('gamma', Parameter) ~ Idx1,
                               !is.na(Idx1) & !grepl('B', Parameter) ~ Idx1),
           AgeIdx = case_when(grepl('^Sigma\\.raw\\[|^xi\\[', Parameter) ~ Idx1,
                              !is.na(Idx2) & !grepl('gamma', Parameter) ~ Idx1,
                              !is.na(Idx2) & grepl('gamma', Parameter) ~ Idx2,
                              is.na(Idx2) & grepl('B', Parameter) ~ Idx1),
           Year = YearIdx + minYear - 1,
           Age = case_when(grepl('^SA\\[|^AD\\[|^s\\.SA\\[|^s\\.AD\\[', Parameter) ~ AgeIdx,
                           grepl('^YAF\\[|^s\\.YAF\\[', Parameter) ~ 1,
                           grepl('^s\\.PY\\[', Parameter) ~ 0,
                           TRUE ~ NA_real_),
           AgeClass = ifelse(grepl('^B\\.age\\[|^B\\.dens\\[|^B\\.veg\\[|^s\\[', Parameter), AgeIdx, NA),
           ParamName = word(Parameter, 1, sep = "\\["))
  
  sum.dat <- sum.dat %>%
    left_join(idx.dat, by = "Parameter")
  
  
  ## Prep parameters for plotting ----------------------------------------------
  
  # set parameter groups for plotting posterior density overlaps
  plot.params <- list(
    CJSbetas = c(paste0('B.age[', 1:nAgeC, ']'),
                 paste0('B.dens[', 1:nAgeC, ']'),
                 paste0('B.veg[', 1:nAgeC, ']')),
    
    CJScovar = c(paste0('Sigma.raw[', 1:nAgeC, ', ', 1:nAgeC, ']')),
                 # expand.grid(t = 1:(ntimes - 1), a = 1:nAgeC) %>%
                 #   mutate(param = paste0('gamma[', t, ', ', a, ']')) %>%
                 #   pull(param),
                 # paste0('xi[', 1:nAgeC, ']')),
    
    CJSobs   = c('mean.p', 'sd.p',
                 paste0('year.p[', 1:ntimes, ']')),
    
    CJSenv   = c(paste0('dens.hat[', 1:(ntimes-1), ']'),
                 paste0('veg.hat[', 1:(ntimes-1), ']')),
    
    CJSsurv  = c(expand.grid(a = 1:nAgeC, t = c(1, 5, 9, 13)) %>% 
                   mutate(param = paste0('s[', a, ', ', t, ']')) %>% 
                   pull(param)),
    
    rates    = c(paste0('b[', 1:(ntimes-1), ']'),
                 paste0('s.PY[', 1:(ntimes-1), ']'))
                 # paste0('s.YAF[', 1:(ntimes-1), ']'),
                 # expand.grid(a = 1:2, t = 1:(ntimes-1)) %>% 
                 #   mutate(param = paste0('s.SA[', a, ', ', t, ']')) %>% 
                 #   pull(param),
                 # expand.grid(a = 1:nAge, t = 1:(ntimes-1)) %>% 
                 #   mutate(param = paste0('s.AD[', a, ', ', t, ']')) %>% 
                 #   pull(param)),
    
    # popsizes = c(paste0('YAF[', 1:ntimes, ']'),
    #              expand.grid(a = 1:2, t = 1:(ntimes-1)) %>% 
    #                mutate(param = paste0('SA[', a, ', ', t, ']')) %>% 
    #                pull(param),
    #              expand.grid(a = 1:nAge, t = 1:(ntimes-1)) %>% 
    #                mutate(param = paste0('AD[', a, ', ', t, ']')) %>% 
    #                pull(param),
    #              paste0('Ntot[', 1:ntimes, ']'))
    )
  
  # set parameters for plotting time series of posterior summaries
  plotTS.rates <- list(
    ParamNames = c('b', 's.PY', 's.YAF', 's.SA', 's.AD'),
    ParamLabels = c('Breeding rate',
                    'Survival to pouch exit',
                    'Survival of young-at-foot (0 yrs)',
                    'Survival of subadults (1-2 yrs)',
                    'Survival of adults (3+ yrs)'))
  
  plotTS.pops <- list(
    ParamNames = c('YAF', 'SA', 'AD', 'Ntot'),
    ParamLabels = c('# of young-at-foot', '# of subadults',
                    '# of adults', 'Total # of females'))
  
  # set plotting colors
  # plot.cols <- paletteer_c("grDevices::Temps", length(modelNames))
  plot.cols <- natparks.pals("Banff", nModels)
  
  
  ## Plot ----------------------------------------------------------------------
  
  # posterior overlaps
  pdf(paste0(plotFolder, "/PosteriorDensities.pdf"), width = 9, height = 6)
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
  pdf(paste0(plotFolder, "/PosteriorSummaries_TSrates.pdf"), width = 8, height = 4)
  for(x in 1:length(plotTS.rates$ParamNames)){
    
    print(
      ggplot(subset(sum.dat, ParamName == plotTS.rates$ParamNames[x] & Year >= minYear & Year <= maxYear), aes(group = Model)) + 
        geom_ribbon(aes(x = Year, ymin = Lower, ymax = Upper, fill = Model), alpha = 1/nModels) + 
        geom_line(aes(x = Year, y = Median, color = Model)) + 
        scale_fill_manual(values = plot.cols) +
        scale_color_manual(values = plot.cols) + 
        scale_x_continuous(breaks = c(minYear:maxYear), labels = c(minYear:maxYear)) + 
        ggtitle(plotTS.rates$ParamLabels[x]) +  
        theme_bw() + theme(panel.grid.minor = element_blank(), 
                           panel.grid.major.y = element_blank(), 
                           axis.text.x = element_text(angle = 45, vjust = 0.5))
    )
    
  }
  dev.off()
  
  # posterior summary time series of population sizes
  pdf(paste0(plotFolder, "/PosteriorSummaries_TSpopsizes.pdf"), width = 8, height = 4)
  for(x in 1:length(plotTS.pops$ParamNames)){
    
    print(
      ggplot(subset(sum.dat, ParamName == plotTS.pops$ParamNames[x] & Year > minYear & Year <= maxYear), aes(group = Model)) + 
        geom_ribbon(aes(x = Year, ymin = Lower, ymax = Upper, fill = Model), alpha = 1/nModels) + 
        geom_line(aes(x = Year, y = Median, color = Model)) + 
        scale_fill_manual(values = plot.cols) +
        scale_color_manual(values = plot.cols) + 
        scale_x_continuous(breaks = c(minYear:maxYear), labels = c(minYear:maxYear)) + 
        ggtitle(plotTS.pops$ParamLabels[x]) +  
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

# test <- compareModels(nAge = 22, nAgeC = 5, ntimes = 17, minYear = 2008,
#                       postPaths = "results/IPM_CJS.rds", modelNames = c("IPM/CJS"),
#                       plotFolder = "figures", returnSumData = TRUE)

