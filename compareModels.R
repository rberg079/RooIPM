#' Compare outputs of different models
#'
#' @param nAgeC 
#' @param ntimes 
#' @param minYear 
#' @param maxYear 
#' @param postPaths 
#' @param postList 
#' @param modelNames 
#' @param plotFolder 
#' @param returnSumData 
#'
#' @returns
#' @export
#'
#' @examples

compareModels <- function(nAgeC, ntimes, minYear, maxYear = minYear,
                          postPaths, postList, modelNames, plotFolder,
                          returnSumData = FALSE){
  
  library(tidyverse)
  library(reshape2)
  library(paletteer)
  
  ## Check models are specified correctly
  if((missing(post.filepaths) & missing(post.list)) |
     (missing(post.filepaths) & missing(post.list))){
    stop("Models have to be specified either via file paths (post.filepaths) or 
         using object names (post.objects).")
  }
  
  ## Make plotting directory if it does not exist already
  if(!dir.exists(plotFolder)){
    dir.create(plotFolder)
  }
  
  ## Count number of models
  nModels <- length(modelNames)
  
  ## Set maxYear if not provided
  if(missing(maxYear)){
    maxYear <- minYear + Tmax - 1
  }
  
  ## Reformat posterior samples
  post.dat <- data.frame()
  for(i in 1:nModels){
    
    # Extract samples for relevant model
    if(!missing(postList)){
      samples <- as.matrix(postList[[i]])
    }else{
      samples <- as.matrix(readRDS(postFilepaths[i]))
    }
    
    # Change format and add to list
    model.dat <- melt(samples)
    colnames(model.dat) <- c("Sample", "Parameter", "Value")
    model.dat$Model <- model.names[i]
    post.dat <- rbind(post.dat, model.dat)
  }
  
  ## Summarize posterior samples into median + 95% CI
  sum.dat <- post.dat %>%
    group_by(Parameter, Model) %>%
    summarise(Median = median(Value, na.rm = TRUE),
              Lower = quantile(Value, probs = 0.025, na.rm = TRUE),
              Upper = quantile(Value, probs = 0.975, na.rm = TRUE),
              .groups = "keep") %>%
    mutate(Parameter = as.character(Parameter)) %>%
    ungroup()
  
  ## Extract and add age and year information
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
           # Age = ifelse(AgeIdx == Amax, paste0(Amax-1, "+"), AgeIdx-1),  # !!!!!! needs attention
           ParamName = word(Parameter, 1, sep = "\\["))
  
  sum.dat <- sum.dat %>%
    left_join(idx.dat, by = "Parameter")
  
  ## Set parameter groups for plotting posterior density overlaps
  plot.params <- list(
    CJSbetas = c(paste0('B.age[', 1:nAgeC, ']'),
                 paste0('B.dens[', 1:nAgeC, ']'),
                 paste0('B.veg[', 1:nAgeC, ']')),
    
    CJScovar = c(expand.grid(t = 1:(ntimes - 1), a = 1:nAgeC) %>%
                   mutate(param = paste0('gamma[', t, ', ', a, ']')) %>%
                   pull(param),
                 paste0('Sigma.raw[', 1:nAgeC, ', ', 1:nAgeC, ']'),
                 paste0('xi[', 1:nAgeC, ']')),
    
    CJSobs   = c('mean.p', 'sd.p',
                 paste0('year.p[', 1:ntimes, ']')),
    
    CJSenv   = c(paste0('dens.hat[', 1:(ntimes-1), ']'),
                 paste0('veg.hat[', 1:(ntimes-1), ']')),
    
    CJSsurv  = c(expand.grid(a = 1:nAgeC, t = 1:(ntimes-1)) %>% 
                   mutate(param = paste0('s[', a, ', ', t, ']')) %>% 
                   pull(param)),
    
    survival    = c(paste0('s.YAF[', 1:(ntimes-1), ']'),
                    expand.grid(a = 1:2, t = 1:(ntimes-1)) %>% 
                      mutate(param = paste0('s.SA[', a, ', ', t, ']')) %>% 
                      pull(param),
                    expand.grid(a = 1:nAge, t = 1:(ntimes-1)) %>% 
                      mutate(param = paste0('s.AD[', a, ', ', t, ']')) %>% 
                      pull(param)),
    
    recruitment = c(paste0('b[', 1:(ntimes-1), ']'),
                    paste0('s.PY[', 1:(ntimes-1), ']')),
    
    popsizes    = c(paste0('YAF[', 1:ntimes, ']'),
                    expand.grid(a = 1:2, t = 1:(ntimes-1)) %>% 
                      mutate(param = paste0('SA[', a, ', ', t, ']')) %>% 
                      pull(param),
                    expand.grid(a = 1:nAge, t = 1:(ntimes-1)) %>% 
                      mutate(param = paste0('AD[', a, ', ', t, ']')) %>% 
                      pull(param),
                    paste0('Ntot[', 1:ntimes, ']')))
  
  ## Set parameters plotting time series of posterior summaries
  plotTS.rates <- list(
    ParamNames = c('b', 's.PY', 's.YAF', 's.SA', 's.AD'),
    ParamLabels = c('Breeding rate',
                    'Survival to pouch exit (PYs)',
                    'Survival of YAFs (0-1 years)',
                    'Survival of SAs (1-2 & 2-3 years)',
                    'Survival of adults (3 years & up)'))
  
  plotTS.pops <- list(
    ParamNames = c('YAF', 'SA', 'AD', 'Ntot'),
    ParamLabels = c('# of YAFs', '# of subadults',
                    '# of adults', 'Total # of females'))
  
  ## Set plotting colors
  plot.cols <- paletteer_c("grDevices::Temps", length(model.names))
  
  ## Plot posterior overlaps
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
  
  ## Plot posterior summary time series for age-specific parameters
  pdf(paste0(plotFolder, "/PosteriorSummariesTS.pdf"), width = 8, height = 8)
  for(x in 1:length(plotTS.pops$ParamNames)){
    
    print(
      ggplot(subset(sum.data, ParamName == plotTS.pops$ParamNames[x] & Year <= maxYear),
             aes(group = Model)) + 
        geom_line(aes(x = Year, y = median, color = Model)) + 
        geom_ribbon(aes(x = Year, ymin = lCI, ymax = uCI, fill = Model), alpha = 1/nModels) + 
        scale_fill_manual(values = plot.cols) +
        scale_color_manual(values = plot.cols) + 
        scale_x_continuous(breaks = c(minYear:maxYear), labels = c(minYear:maxYear)) + 
        facet_wrap(~ Age, ncol = 1, scales = "free_y") + 
        ggtitle(plotTS.pops$ParamLabels[x]) +  
        theme_bw() + theme(panel.grid.minor = element_blank(), 
                           panel.grid.major.y = element_blank(), 
                           axis.text.x = element_text(angle = 45, vjust = 0.5))
    )
  }
  dev.off()
  
  ## Optional: return summary data
  if(returnSumData){
    return(sum.data)
  }
}

test <- compareModels(nAgeC = 5, ntimes = 17, minYear = 2008,
                      postFilepaths = "results/IPM_CJS.rds", modelNames = c("IPM/CJS"),
                      plotFolder = "figures", returnSumData = TRUE)
test

