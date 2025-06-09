

plotLTRE_random <- function(LTREresults, nAge = 19)
  
# for testing purposes
LTREresults <- readRDS('results/LTREresults_random.rds')
# plotFolder = c("figures")
nAge = 19


## Set up --------------------------------------------------------------------

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
  dplyr::filter(variable %in% c("Bt", "Ra_sum",
                                "sYAF", "sSA", "sAD_sum",
                                "pYAF", "pSA", "pAD_sum")) %>% 
  dplyr::mutate(type = dplyr::case_when(variable == "Bt" ~ "Birth rate",
                                        variable == "Ra_sum" ~ "Survival of jellybeans",
                                        variable == "sYAF" ~ "Survival of young-at-foot",
                                        variable == "sSA" ~ "Survival of subadults",
                                        variable == "sAD_sum" ~ "Survival of adults",
                                        variable == "pYAF" ~ "Proportion of young-at-foot",
                                        variable == "pSA" ~ "Proportion of subadults",
                                        variable == "pAD_sum" ~ "Proportion of adults"))

# make ordered list of parameter types
typeList <- c("Birth rate", "Survival of jellybeans",
              "Survival of young-at-foot", "Survival of subadults", "Survival of adults",
              "Proportion of young-at-foot", "Proportion of subadults", "Proportion of adults")

# order factor levels
contData_sum$type <- factor(contData_sum$type, levels = typeList)


## Plot contributions (violin plots) -----------------------------------------

# plot colours
temp.colours <- paletteer::paletteer_c("grDevices::Temps", length(unique(contData_sum$type)))

# plot colours
temp.colours <- paletteer::paletteer_c("grDevices::Temps", length(unique(sum.data$type)))
# plot.colours <- c("#047993FF", "#005F94FF", temp.colours[1:3], rep(temp.colours[4], 2), temp.colours[5:6])
plot.colours <- temp.colours

# summed estimates for all parameters
addline_format <- function(x,...){
  gsub('\\s','\n',x)
}

p.sum <- contData_sum %>% 
  dplyr::filter(Contribution < 0.5) %>% 
  ggplot(aes(x = type, y = Contribution, group = type)) +
  geom_violin(aes(fill = type, colour = type), alpha = 0.5, scale = "width", draw_quantiles = 0.5) +
  geom_hline(yintercept = 0, colour = "grey70", linetype = "dashed") +
  ylab("Contribution") +
  xlab("") +
  scale_fill_manual(values = plot.colours) +
  scale_colour_manual(values = plot.colours) +
  scale_x_discrete(labels = addline_format(typeList)) + # WITH ADDLINE!
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.title = element_text(size = 10))

# reproductive success panel
R.colours <- c(plot.colours[1], rep(plot.colours[2], 18))
names(R.colours) <- c("Bt", paste0("Ra_", 2:nAge))

p.R <- ggplot(subset(contData, variable %in% c("Bt", paste0("Ra_", 2:nAge))),
                     aes(x = variable, y = Contribution, group = variable)) +
  geom_violin(fill = R.colours, colour = R.colours, alpha = 0.5, scale = "width", draw_quantiles = 0.5) +
  geom_hline(yintercept = 0, colour = "grey70", linetype = "dashed") +
  ylab("Contribution") +
  xlab("") +
  scale_x_discrete(labels = expression(B, R[2], R[3], R[4], R[5], R[6], R[7], R[8], R[9], R[10],
                                       R[11], R[12], R[13], R[14], R[15], R[16], R[17], R[18], R[19])) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.title = element_text(size = 10))

# survival panel
S.colours <- c(plot.colours[3:4], rep(plot.colours[5], 18))
names(S.colours) <- c("sYAF", "sSA", paste0("sAD_", 2:nAge))

p.S <- ggplot(subset(contData, ))








#-----------------------------------#
# Plot contributions - Violin plots #
#-----------------------------------#

if(!HazardRates){
  ## Survival panel
  p.S <- ggplot(subset(contData, variable2 %in% paste0("S_", 1:Amax)), aes(x = variable2, y = Contribution, fill = Season)) + 
    geom_violin(color = plot.colors[1], alpha = 0.5, scale = 'width', draw_quantiles = 0.5) + 
    geom_hline(yintercept = 0, color = "grey70", linetype = "dashed") + 
    ylab("Contribution") + 
    xlab('') + 
    scale_x_discrete(labels = expression(S[1], S[2], S[3], S[4], S[5])) +
    scale_fill_manual(values = c("white", plot.colors[1])) +
    theme_bw() + 
    theme(legend.position = c(0.8, 0.8), 
          legend.key.size = unit(0.5, 'cm'),
          legend.key.height = unit(0.5, 'cm'),
          legend.key.width = unit(0.5, 'cm'), 
          legend.title = element_blank(), 
          legend.text = element_text(size = 8),
          panel.grid = element_blank(), 
          axis.text.x = element_text(size = 12), 
          axis.title = element_text(size = 12))
}else{
  ## Harvest mortality panel
  p.mH <- ggplot(subset(contData, variable2 %in% paste0("mH_", 1:Amax)), aes(x = variable2, y = Contribution, fill = Season)) + 
    geom_violin(color = plot.colors[1], alpha = 0.5, scale = 'width', draw_quantiles = 0.5) + 
    geom_hline(yintercept = 0, color = "grey70", linetype = "dashed") + 
    ylab("Contribution") + 
    xlab('') + 
    scale_x_discrete(labels = expression(m[1]^H, m[2]^H, m[3]^H, m[4]^H, m[5]^H)) +
    scale_fill_manual(values = c("white", plot.colors[1])) +
    theme_bw() + 
    theme(legend.position = c(0.8, 0.8), 
          legend.key.size = unit(0.5, 'cm'),
          legend.key.height = unit(0.5, 'cm'),
          legend.key.width = unit(0.5, 'cm'), 
          legend.title = element_blank(), 
          legend.text = element_text(size = 8),
          panel.grid = element_blank(), 
          axis.text.x = element_text(size = 12), 
          axis.title = element_text(size = 12))
  
  ## Natural mortality panel
  p.mO <- ggplot(subset(contData, variable %in% paste0("mO_", 1:Amax)), aes(x = variable, y = Contribution, group = variable)) + 
    geom_violin(fill = plot.colors[2], color = plot.colors[2], alpha = 0.5, scale = 'width', draw_quantiles = 0.5) + 
    geom_hline(yintercept = 0, color = "grey70", linetype = "dashed") + 
    ylab("Contribution") + 
    xlab('') + 
    scale_x_discrete(labels = expression(m[1]^O, m[2]^O, m[3]^O, m[4]^O, m[5]^O)) + 
    theme_bw() + 
    theme(legend.position = 'none', panel.grid = element_blank(), axis.text.x = element_text(size = 12), axis.title = element_text(size = 12))
}

## Population structure panel
if(PopStructure){
  p.n <- ggplot(subset(contData, variable %in% paste0("n_", 1:Amax)), aes(x = variable, y = Contribution, group = variable)) + 
    geom_violin(fill = plot.colors[length(plot.colors)], color = plot.colors[length(plot.colors)], alpha = 0.5, scale = 'width', draw_quantiles = 0.5) + 
    geom_hline(yintercept = 0, color = "grey70", linetype = "dashed") + 
    ylab("Contribution") + 
    xlab('') + 
    scale_x_discrete(labels = expression(n[1], n[2], n[3], n[4], n[5])) + 
    theme_bw() + 
    theme(legend.position = 'none', panel.grid = element_blank(), axis.text.x = element_text(size = 12), axis.title = element_text(size = 12))
}else{
  p.n <- ggplot(subset(contData, variable %in% paste0("N_", 1:Amax)), aes(x = variable, y = Contribution, group = variable)) + 
    geom_violin(fill = plot.colors[length(plot.colors)], color = plot.colors[length(plot.colors)], alpha = 0.5, scale = 'width', draw_quantiles = 0.5) + 
    geom_hline(yintercept = 0, color = "grey70", linetype = "dashed") + 
    ylab("Contribution") + 
    xlab('') + 
    scale_x_discrete(labels = expression(N[1], N[2], N[3], N[4], N[5])) + 
    theme_bw() + 
    theme(legend.position = 'none', panel.grid = element_blank(), axis.text.x = element_text(size = 12), axis.title = element_text(size = 12))
}


## Combine panels and save to pdf
pdf(paste0("Plots/RedFoxIPM_randomLTRE_", ifelse(HazardRates, "MHR", "SP"), "_sum.pdf"), width = 10, height = 6)
print(
  p.sum
)
dev.off()

if(HazardRates){
  
  pdf(paste0("Plots/RedFoxIPM_randomLTRE_", ifelse(HazardRates, "MHR", "SP"), "_age.pdf"), width = 7, height = 8)
  print(
    (p.mH + labs(tag = 'a)') | p.mO + labs(tag = 'b)')) / (p.Psi + labs(tag = 'c)')| p.rho + labs(tag = 'd)')) / (p.n  + labs(tag = 'e)') | plot_spacer())
  )
  dev.off()
  
}else{
  pdf(paste0("Plots/RedFoxIPM_randomLTRE_", ifelse(HazardRates, "MHR", "SP"), "_age.pdf"), width = 7, height = 6)
  print(
    (p.S + labs(tag = 'a)')| p.Psi + labs(tag = 'b)')) / (p.rho  + labs(tag = 'c)')| p.n + labs(tag = 'd)'))
  )
  dev.off()
}

## Return list of plots
plotList <- c(paste0("Plots/RedFoxIPM_randomLTRE_", ifelse(HazardRates, "MHR", "SP"), "_sum.pdf"),
              paste0("Plots/RedFoxIPM_randomLTRE_", ifelse(HazardRates, "MHR", "SP"), "_age.pdf"))

return(plotList)