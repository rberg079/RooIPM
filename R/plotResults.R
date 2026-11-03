# 20 April 2026
# Plot results of IPM

## Set up ----------------------------------------------------------------------

library(tidyverse)
library(data.table)
library(paletteer)
library(patchwork)
library(scales)

nYear <- 18
nAge  <- 19

# load results
out.mcmc <- readRDS('results/IPM_CJSen_RSen_AB_DynDens_dCJS_Dave12_S33_R22_B22_muFunE_centerAge_areaHR.rds')
out.mat <- do.call(rbind, lapply(out.mcmc, as.matrix))


## Population size -------------------------------------------------------------

# parameters to include
table.params <- c(
  paste0('nYF[', 1:nYear, ']'),
  paste0('nSA[', 1:nYear, ']'),
  paste0('nAD[', rep(1:nAge, times = nYear), ', ', rep(1:nYear, each = nAge), ']'),
  paste0('nTOT[', 1:nYear, ']'))

# build summary dataframe
post.table <- data.frame(Parameter = table.params, Estimate = NA)

for(i in 1:length(table.params)){
  est <- out.mat[, table.params[i]]
  post.table$Estimate[i] <- paste0(round(median(est, na.rm = T), digits = 2), ' [',
                                   round(quantile(est, 0.025, na.rm = T), digits = 2), ', ',
                                   round(quantile(est, 0.975, na.rm = T), digits = 2), ']')
}

# indices
nYF <- grep("^nYF\\[", colnames(out.mcmc[[1]])); nYF
nSA <- grep("^nSA\\[", colnames(out.mcmc[[1]])); nSA
nAD <- grep("^nAD\\[", colnames(out.mcmc[[1]])); nAD
nTOT <- grep("^nTOT\\[", colnames(out.mcmc[[1]])); nTOT

# select variable!
var <- nTOT

popDF <- data.frame(
  Year = 1:length(var),
  Mean = apply(out.mat[, var, drop = FALSE], 2, mean, na.rm = TRUE),
  Lower = apply(out.mat[, var, drop = FALSE], 2, quantile, probs = 0.025, na.rm = TRUE),
  Upper = apply(out.mat[, var, drop = FALSE], 2, quantile, probs = 0.975, na.rm = TRUE)
)

# population plot
pop <- popDF %>% 
  mutate(Year = Year + 2007) %>% 
  ggplot(aes(x = Year, y = Mean)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "#62414B", alpha = 0.2) +
  geom_line(color = "#62414B", linewidth = 0.8) +
  scale_x_continuous(limits = c(2008, 2025),
                     breaks = c(2008, 2010, 2012, 2014, 2016, 2018, 2020, 2022, 2024)) +
  scale_y_continuous(limits = c(100, 450),
                     breaks = c(100, 200, 300, 400, 500)) +
  labs(y = "Population size") +
  theme_bw() +
  theme(
        axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none"
        ); pop

# ggsave("figures/results/final!/nTOT.png", width = 18.0, height = 10.0, units = c("cm"), dpi = 600)


## Population size (by group) --------------------------------------------------

# parameters to include
nYF_idx  <- grep("^nYF\\[", colnames(out.mat))
nSA_idx  <- grep("^nSA\\[", colnames(out.mat))
nTOT_idx <- grep("^nTOT\\[", colnames(out.mat))
nYear    <- length(nYF_idx)

mat_2to9 <- matrix(NA, nrow = nrow(out.mat), ncol = nYear)
mat_10up <- matrix(NA, nrow = nrow(out.mat), ncol = nYear)

for(y in 1:nYear) {
  cols_all  <- grep(paste0("^nAD\\[[0-9]+,\\s*", y, "\\]$"), colnames(out.mat)) # all adults
  cols_2to9 <- grep(paste0("^nAD\\[[2-9],\\s*", y, "\\]$"), colnames(out.mat))  # 2-9 years
  cols_10up <- setdiff(cols_all, cols_2to9)                                     # 10+ years
  
  # sum across adult age groups
  mat_2to9[, y] <- rowSums(out.mat[, cols_2to9, drop = FALSE], na.rm = TRUE)
  mat_10up[, y] <- rowSums(out.mat[, cols_10up, drop = FALSE], na.rm = TRUE)
}

# build summary dataframes
df_YF   <- data.frame(Year = 1:nYear, Group = "0 years",
                      Mean  = apply(out.mat[, nYF_idx, drop = FALSE], 2, mean, na.rm = TRUE),
                      Lower = apply(out.mat[, nYF_idx, drop = FALSE], 2, quantile, probs = 0.025, na.rm = TRUE),
                      Upper = apply(out.mat[, nYF_idx, drop = FALSE], 2, quantile, probs = 0.975, na.rm = TRUE))

df_SA   <- data.frame(Year = 1:nYear, Group = "1 year",
                      Mean  = apply(out.mat[, nSA_idx, drop = FALSE], 2, mean, na.rm = TRUE),
                      Lower = apply(out.mat[, nSA_idx, drop = FALSE], 2, quantile, probs = 0.025, na.rm = TRUE),
                      Upper = apply(out.mat[, nSA_idx, drop = FALSE], 2, quantile, probs = 0.975, na.rm = TRUE))

df_2to9 <- data.frame(Year = 1:nYear, Group = "2-9 years",
                      Mean  = apply(mat_2to9, 2, mean, na.rm = TRUE),
                      Lower = apply(mat_2to9, 2, quantile, probs = 0.025, na.rm = TRUE),
                      Upper = apply(mat_2to9, 2, quantile, probs = 0.975, na.rm = TRUE))

df_10up <- data.frame(Year = 1:nYear, Group = "10+ years",
                      Mean  = apply(mat_10up, 2, mean, na.rm = TRUE),
                      Lower = apply(mat_10up, 2, quantile, probs = 0.025, na.rm = TRUE),
                      Upper = apply(mat_10up, 2, quantile, probs = 0.975, na.rm = TRUE))

df_TOT  <- data.frame(Year = 1:nYear, Group = "Total",
                      Mean  = apply(out.mat[, nTOT_idx, drop = FALSE], 2, mean, na.rm = TRUE),
                      Lower = apply(out.mat[, nTOT_idx, drop = FALSE], 2, quantile, probs = 0.025, na.rm = TRUE),
                      Upper = apply(out.mat[, nTOT_idx, drop = FALSE], 2, quantile, probs = 0.975, na.rm = TRUE))

# combine everything
df <- rbind(df_YF, df_SA, df_2to9, df_10up, df_TOT)

# # pick colours
# cols <- c(
#   "Young-at-foot (0 years)" = "#8E5EA2",
#   "Subadults (1 year)"      = "#277DA1",
#   "Adults (2-9 years)"      = "#F9C74F",
#   "Adults (10+ years)"      = "#D62828",
#   "Total female population" = "#62414B"
# )

cols <- c(
  "0 years"   = "#AC2020",
  "1 year"    = "#DD962C",
  "2-9 years" = "#E3CA78",
  "10+ years" = "#6F9A9B",
  "Total"     = "black"
)

# population plot
popA <- df %>% 
  mutate(Year = Year + 2007,
         Group = factor(Group,
                        levels = c("0 years",
                                   "1 year",
                                   "2-9 years",
                                   "10+ years",
                                   "Total"))) %>% 
  ggplot(aes(x = Year, y = Mean, group = Group, colour = Group, fill = Group)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.2, colour = NA) + 
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = cols) +
  scale_fill_manual(values = cols) +
  scale_x_continuous(limits = c(2008, 2025),
                     breaks = c(2008, 2010, 2012, 2014, 2016, 2018, 2020, 2022, 2024)) +
  scale_y_continuous(limits = c(0, 450), breaks = scales::pretty_breaks()) +
  labs(y = "Population size", colour = "Age group", fill = "Age group") +
  theme_bw(); popA

# ggsave("figures/results/final!/allNs3.png", width = 18.0, height = 10.0, units = c("cm"), dpi = 600)


## Survival --------------------------------------------------------------------

# indices
S_idx  <- grep("^S\\[", colnames(out.mat))

# build summary dataframe
df <- expand.grid(Age = 1:13, Year = 1:17)
df$Mean <- apply(out.mat[, S_idx, drop = FALSE], 2, mean, na.rm = TRUE)
df$Lower <- apply(out.mat[, S_idx, drop = FALSE], 2, quantile, probs = 0.025, na.rm = TRUE)
df$Upper <- apply(out.mat[, S_idx, drop = FALSE], 2, quantile, probs = 0.975, na.rm = TRUE)

# pick colours
cols <- c(
  "0"   = "#332288",
  "1"   = "#0072B2",
  "2"   = "#009E73",
  "4"   = "#F0E442",
  "6"   = "#E69F00",
  "8"   = "#D55E00",
  "10"  = "#CC3311",
  "12+" = "#CC79A7"
)

# pick linetypes
lts <- c(
  "0"   = "solid",
  "1"   = "dashed",
  "2"   = "solid",
  "4"   = "dashed",
  "6"   = "solid",
  "8"   = "dashed",
  "10"  = "solid",
  "12+" = "dashed"
)

# plot
surv <- df %>%
  mutate(Year = Year + 2007,
         Age  = factor(Age-1,
                       levels = c(0, 1, 2, 4, 6, 8, 10, 12),
                       labels = c("0", "1", "2", "4", "6", "8", "10", "12+"))) %>%
  filter(!is.na(Age)) %>% 
  ggplot(aes(x = Year, y = Mean, group = Age, colour = Age, linetype = Age)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = Age), alpha = 0.2, colour = NA) +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = cols) +
  scale_fill_manual(values = cols) +
  scale_linetype_manual(values = lts) +
  scale_x_continuous(limits = c(2008, 2024),
                     breaks = c(2008, 2010, 2012, 2014, 2016, 2018, 2020, 2022, 2024)) +
  scale_y_continuous(limits = c(0, 1), breaks = pretty_breaks()) +
  labs(x = "Year", y = "Survival", colour = "Age", fill = "Age") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank()
    ); surv

# ggsave("figures/results/final!/survival.png", width = 18.0, height = 10.0, units = c("cm"), dpi = 600)


## Reproductive success --------------------------------------------------------

# indices
BR_idx  <- grep("^BR\\[", colnames(out.mat))
sPY_idx <- grep("^sPY\\[", colnames(out.mat))

# extract matrices
BR  <- out.mat[, BR_idx, drop = FALSE]
sPY <- out.mat[, sPY_idx, drop = FALSE]

# compute reproductive output
R <- 0.5 * BR * sPY

# build summary dataframe
df <- expand.grid(Age = 1:19, Year = 1:17)
df$Mean <- apply(R, 2, mean, na.rm = TRUE)
df$Lower <- apply(R, 2, quantile, probs = 0.025, na.rm = TRUE)
df$Upper <- apply(R, 2, quantile, probs = 0.975, na.rm = TRUE)

# # build summary dataframe for BIRTH RATE
# df <- expand.grid(Age = 1:19, Year = 1:17)
# df$Mean <- apply(out.mat[, BR_idx, drop = FALSE], 2, mean, na.rm = TRUE)
# df$Lower <- apply(out.mat[, BR_idx, drop = FALSE], 2, quantile, probs = 0.025, na.rm = TRUE)
# df$Upper <- apply(out.mat[, BR_idx, drop = FALSE], 2, quantile, probs = 0.975, na.rm = TRUE)
 
# # build summary dataframe for PY SURVIVAL
# df <- expand.grid(Age = 1:19, Year = 1:17)
# df$Mean <- apply(out.mat[, sPY_idx, drop = FALSE], 2, mean, na.rm = TRUE)
# df$Lower <- apply(out.mat[, sPY_idx, drop = FALSE], 2, quantile, probs = 0.025, na.rm = TRUE)
# df$Upper <- apply(out.mat[, sPY_idx, drop = FALSE], 2, quantile, probs = 0.975, na.rm = TRUE)

# pick colours
cols <- c(
  "2"   = "#009E73",
  "4"   = "#F0E442",
  "6"   = "#E69F00",
  "8"   = "#D55E00",
  "10"  = "#CC3311",
  "12+" = "#CC79A7"
)

# pick linetypes
lts <- c(
  "2"   = "solid",
  "4"   = "dashed",
  "6"   = "solid",
  "8"   = "dashed",
  "10"  = "solid",
  "12+" = "dashed"
)

# plot
rs <- df %>%
  mutate(Year = Year + 2007,
         Age  = factor(Age,
                       levels = c(2, 4, 6, 8, 10, 12),
                       labels = c("2", "4", "6", "8", "10", "12+"))) %>%
  filter(!is.na(Age)) %>% 
  ggplot(aes(x = Year, y = Mean, group = Age, colour = Age, linetype = Age)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = Age), alpha = 0.2, colour = NA, show.legend = F) +
  geom_line(linewidth = 0.8, show.legend = F) + # position = position_dodge(width = 0.2), 
  scale_colour_manual(values = cols) +
  scale_fill_manual(values = cols) +
  scale_linetype_manual(values = lts) +
  scale_x_continuous(breaks = c(2008, 2010, 2012, 2014, 2016, 2018, 2020, 2022, 2024)) +
  scale_y_continuous(limits = c(0, 1), breaks = pretty_breaks()) +
  labs(x = "Year", y = "Reproductive success", colour = "Age", fill = "Age") +
  # labs(x = "Year", y = "Birth rate", colour = "Age", fill = "Age") +
  # labs(x = "Year", y = "Survival of pouch young", colour = "Age", fill = "Age") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_blank()
  ); rs

# ggsave("figures/results/final!/RS.png", width = 18.0, height = 10.0, units = c("cm"), dpi = 600)

# # combine birth rate & PY survival
# (br / rs) +
#   plot_layout(guides = "collect") +
#   theme(legend.position = "right")

# ggsave("figures/results/final!/birth&PYsurv.png", width = 18.0, height = 18.0, units = c("cm"), dpi = 600)

# # summaries to report
# R_array <- array(R, dim = c(nrow(R), 19, 17))
# R_age <- apply(R_array, c(1, 2), mean, na.rm = TRUE)
# 
# df_age <- data.frame(Age   = 1:19,
#   Mean  = apply(R_age, 2, mean),
#   Lower = apply(R_age, 2, quantile, 0.025),
#   Upper = apply(R_age, 2, quantile, 0.975))
# 
# R_year <- apply(R_array, c(1, 3), mean, na.rm = TRUE)
# 
# df_year <- data.frame(Year  = 1:17,
#   Mean  = apply(R_year, 2, mean),
#   Lower = apply(R_year, 2, quantile, 0.025),
#   Upper = apply(R_year, 2, quantile, 0.975))


## Covariate effects -----------------------------------------------------------

inv_logit <- function(x) 1 / (1 + exp(-x))

keep_ages <- c(1, 11, 12, 13)
age_labels <- c("0", "10", "11", "12+")

# extract betas
bDSy <- out.mat[, grep("^BetaD\\.Sy", colnames(out.mat))]
bDSp <- out.mat[, grep("^BetaD\\.Sp", colnames(out.mat))]
bDSo <- out.mat[, grep("^BetaD\\.So", colnames(out.mat))]
bVSy <- out.mat[, grep("^BetaV\\.Sy", colnames(out.mat))]
bVSp <- out.mat[, grep("^BetaV\\.Sp", colnames(out.mat))]
bVSo <- out.mat[, grep("^BetaV\\.So", colnames(out.mat))]

bDR <- out.mat[, grep("^BetaD\\.R", colnames(out.mat))]

# extract baseline intercepts
muS  <- paste0("Mu.S[", keep_ages, "]")
base <- out.mat[, muS, drop = FALSE]

muR <- grep("^Mu\\.R\\[", colnames(out.mat))
baseR <- rowMeans(out.mat[, muR, drop = FALSE])

# covariate sequence
x <- seq(-2, 2, length.out = 50)

# build summary dfs
make_df <- function(betas_young, betas_prime, betas_old, covariate){
  results <- list()
  
  for(i in seq_along(keep_ages)){
    age <- keep_ages[i]
    label <- age_labels[i]
    intercepts <- base[, i]
    
    # assign the right beta
    if(age == 1){
      betas <- betas_young
    }else if(between(age, 2, 9)){
      betas <- betas_prime
    }else{
      betas <- betas_old
    }
    
    for(xx in x){
      preds <- inv_logit(intercepts + betas * xx)
      results[[length(results) + 1]] <- data.frame(
        Age       = label,
        x         = xx,
        Mean      = mean(preds),
        Lower     = quantile(preds, 0.025),
        Upper     = quantile(preds, 0.975),
        covariate = covariate
      )
    }
  }
  do.call(rbind, results)
}

# build summary df for PY
make_py_df <- function(betas, covariate){
  results <- list()
  for(xx in x){
    preds <- inv_logit(baseR + betas * xx)
    results[[length(results) + 1]] <- data.frame(
      Age       = "PY",
      x         = xx,
      Mean      = mean(preds),
      Lower     = quantile(preds, 0.025),
      Upper     = quantile(preds, 0.975),
      covariate = covariate
    )
  }
  do.call(rbind, results)
}

df <- rbind(
  make_df(bDSy, bDSp, bDSo, "Population density"),
  make_df(bVSy, bVSp, bVSo, "Forage availability"),
  make_py_df(bDR, "Population density")
)

# pick colours
cols <- c(
  "PY"  = "#A578BA",
  "0"   = "#332288",
  "10"  = "#CC3311",
  "11"  = "#696969",
  "12+" = "#CC79A7"
)

# plot
df %>%
  mutate(Age = factor(Age, levels = c("PY", "0", "10", "11", "12+")),
         covariate = factor(covariate, levels = c("Population density", "Forage availability")),
         line = ifelse(Age %in% c("PY", "0") & covariate == "Population density", "solid", "dotted")) %>%
  ggplot(aes(x = x, y = Mean, colour = Age)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = Age), alpha = 0.2, colour = NA) +
  geom_line(aes(linetype = line), linewidth = 0.8) +
  facet_wrap(~ covariate, scales = "free_x") +
  scale_colour_manual(values = cols) +
  scale_fill_manual(values = cols) +
  scale_linetype_identity() +
  scale_y_continuous(breaks = pretty_breaks()) +
  labs(x = "Scaled covariate value", y = "Survival",
       colour = "Age", fill = "Age") +
  guides(linetype = "none") +
  theme_bw()

# ggsave("figures/results/final!/coveffects.png", width = 18.0, height = 10.0, units = c("cm"), dpi = 600)


## Covariate values ------------------------------------------------------------

# indices
D_idx  <- grep("^dens\\.true", colnames(out.mat))#[1:17]
V_idx  <- grep("^veg\\.true", colnames(out.mat))

# extract means and 95% credible intervals
dens_mean <- apply(out.mat[, D_idx, drop = FALSE], 2, mean, na.rm = TRUE)
dens_lci  <- apply(out.mat[, D_idx, drop = FALSE], 2, quantile, probs = 0.025, na.rm = TRUE)
dens_uci  <- apply(out.mat[, D_idx, drop = FALSE], 2, quantile, probs = 0.975, na.rm = TRUE)

veg_mean  <- apply(out.mat[, V_idx, drop = FALSE], 2, mean, na.rm = TRUE)
veg_lci   <- apply(out.mat[, V_idx, drop = FALSE], 2, quantile, probs = 0.025, na.rm = TRUE)
veg_uci   <- apply(out.mat[, V_idx, drop = FALSE], 2, quantile, probs = 0.975, na.rm = TRUE)

# scale the density mean, then apply the exact same center and scale to the CIs
dens_scaled <- scale(dens_mean)
center_val  <- attr(dens_scaled, "scaled:center")
scale_val   <- attr(dens_scaled, "scaled:scale")

dens_lci_scaled <- (dens_lci - center_val) / scale_val
dens_uci_scaled <- (dens_uci - center_val) / scale_val

# build summary dataframes
df_dens <- data.frame(
  Year = 1:18, 
  value = dens_scaled[, 1], 
  lci = dens_lci_scaled, 
  uci = dens_uci_scaled, 
  covariate = "Population density"
)

df_veg <- data.frame(
  Year = 1:17, 
  value = veg_mean, 
  lci = veg_lci, 
  uci = veg_uci, 
  covariate = "Forage availability"
)

df <- rbind(df_dens, df_veg)

# pick colours
line_cols <- c(
  "Population density"  = "#5B5365",
  "Forage availability" = "#7D9570"
)

fill_cols <- c(
  "Population density"  = adjustcolor("#5B5365", alpha.f = 0.10),
  "Forage availability" = adjustcolor("#7D9570", alpha.f = 0.20)
)

# plot
covs <- df %>%
  mutate(Year = Year + 2007,
         covariate = factor(covariate, levels = c("Population density", "Forage availability"))) %>%
  ggplot(aes(x = Year, y = value, colour = covariate, fill = covariate)) +
  geom_hline(yintercept = 0, colour = "grey40") +
  geom_ribbon(aes(ymin = lci, ymax = uci), colour = NA) + 
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = line_cols) +
  scale_fill_manual(values = fill_cols) +
  scale_x_continuous(limits = c(2008, 2024),
                     breaks = c(2008, 2010, 2012, 2014, 2016, 2018, 2020, 2022, 2024)) +
  labs(y = "Scaled covariate value", colour = "Covariate", fill = "Covariate") +
  theme_bw(); covs

# ggsave("figures/results/final!/covsVStime.png", width = 18.0, height = 10.0, units = c("cm"), dpi = 600)

# combine with survival & reproductive output
(surv / rs / covs) +
  plot_layout(guides = "collect") +
  theme(legend.position = "right")

# ggsave("figures/results/final!/surv&rs&covs.png", width = 20.0, height = 22.0, units = c("cm"), dpi = 600)

# plot observed vs true density
myData$densE[1] <- NA

D_idx_all <- grep("^dens\\.true", colnames(out.mat))[1:18]

true_dens_mean <- apply(out.mat[, D_idx_all, drop = FALSE], 2, mean, na.rm = TRUE)
true_dens_lci  <- apply(out.mat[, D_idx_all, drop = FALSE], 2, quantile, probs = 0.025, na.rm = TRUE)
true_dens_uci  <- apply(out.mat[, D_idx_all, drop = FALSE], 2, quantile, probs = 0.975, na.rm = TRUE)

# combine & calculate 95% CrIs for observed density 
# (densE is standard error: Mean +/- 1.96 * SE)
df_obs <- data.frame(
  Year = (1:18) + 2007,
  True_Mean = true_dens_mean,
  True_LCI  = true_dens_lci,
  True_UCI  = true_dens_uci,
  Obs_Mean  = myData$dens,
  Obs_SE    = myData$densE
) %>%
  mutate(
    Obs_LCI = Obs_Mean - (1.96 * Obs_SE),
    Obs_UCI = Obs_Mean + (1.96 * Obs_SE)
  )

# plot
obs_vs_true <- df_obs %>%
  filter(Year > 2008) %>%
  ggplot(aes(x = Year)) +
  
  geom_ribbon(aes(ymin = Obs_LCI, ymax = Obs_UCI, fill = "Observed (Dave's)"), alpha = 0.15) +
  geom_line(aes(y = Obs_Mean, colour = "Observed (Dave's)"), linewidth = 0.8, alpha = 0.8) +
  # geom_errorbar(aes(ymin = Obs_LCI, ymax = Obs_UCI, colour = "Observed"), width = 0.3) +
  # geom_point(aes(y = Obs_Mean, colour = "Observed"), size = 2) +
  
  geom_ribbon(aes(ymin = True_LCI, ymax = True_UCI, fill = "Estimated (IPM's)"), alpha = 0.2) +
  geom_line(aes(y = True_Mean, colour = "Estimated (IPM's)"), linewidth = 0.8) +
  scale_colour_manual(values = c("Estimated (IPM's)" = "#335B5B", "Observed (Dave's)" = "#D55E00")) +
  scale_fill_manual(values = c("Estimated (IPM's)" = "#335B5B", "Observed (Dave's)" = "#D55E00")) +
  scale_x_continuous(breaks = seq(2008, 2024, by = 2)) +
  labs(x = "Year", y = "Population density", colour = "Source", fill = "Source") +
  theme_bw(); obs_vs_true

# ggsave("figures/results/final!/obsVStrueDens.png", width = 18.0, height = 10.0, units = c("cm"), dpi = 600)


## Lambda ----------------------------------------------------------------------

# load results
paramSamples <- readRDS('results/paramSamples.rds')

# extract lambda
lambda_mat <- paramSamples$t$lambda

# summary stats
mean <- apply(lambda_mat, 1, mean)
quantile(mean, c(0.025, 0.5, 0.975))

# build summary dataframe
lamDF <- data.frame(Year  = 1:17,
                    Mean  = apply(lambda_mat, 2, mean),
                    Lower = apply(lambda_mat, 2, quantile, 0.025),
                    Upper = apply(lambda_mat, 2, quantile, 0.975))

lambda <- lamDF %>%
  mutate(Year = Year + 2007) %>%
  ggplot(aes(x = Year, y = Mean)) +
  geom_hline(yintercept = 1, colour = "grey60") +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.2, fill = "#8D6B48") +
  geom_line(colour = "#8D6B48", linewidth = 0.8) +
  scale_x_continuous(limits = c(2008, 2025),
                     breaks = c(2008, 2010, 2012, 2014, 2016, 2018, 2020, 2022, 2024)) +
  # scale_x_continuous(limits = c(2008, 2025),
  #                    breaks = c(2009, 2011, 2013, 2015, 2017, 2019, 2021, 2023, 2025)) +
  scale_y_continuous(limits = c(NA, 1.31),
                     breaks = c(0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3)) +
  labs(x = "Year", y = expression("Population growth rate" ~ (lambda))) +
  theme_bw(); lambda

# ggsave("figures/results/final!/lambda.png", width = 18.0, height = 10.0, units = c("cm"), dpi = 600)

# combine with pop size
pop / lambda
# pAGE / lambda

# overlay with pop size
pop_min <- 100
pop_max <- 450
lam_min <- 0.6
lam_max <- 1.4
s_fac   <- (pop_max - pop_min) / (lam_max - lam_min)
s_shift <- pop_min - (lam_min * s_fac)

popDF <- popDF %>% mutate(Year = Year + 2007)
lamDF <- data.frame(Year  = 1:17 + 2007,
                    Mean  = apply(lambda_mat, 2, mean) * s_fac + s_shift,
                    Lower = apply(lambda_mat, 2, quantile, 0.025) * s_fac + s_shift,
                    Upper = apply(lambda_mat, 2, quantile, 0.975) * s_fac + s_shift)

ggplot() +
  geom_hline(yintercept = 1 * s_fac + s_shift, colour = "grey60", linetype = "dashed") +
  # geom_hline(yintercept = 1.006 * s_fac + s_shift, colour = "#AC2020", linewidth = 0.6) +
  # annotate("text", x = 2023, y = 1.006 * s_fac + s_shift, label = expression("Mean" ~ lambda ~ "(1.006)"),
  #          hjust = 0, vjust = -0.6, size = 3.2, colour = "#AC2020") +
  geom_ribbon(data = popDF, aes(x = Year, ymin = Lower, ymax = Upper), fill = "black", alpha = 0.2) +
  geom_line(data = popDF, aes(x = Year, y = Mean), color = "black", linewidth = 0.8) +
  geom_ribbon(data = lamDF, aes(x = Year, ymin = Lower, ymax = Upper), alpha = 0.2, fill = "#7A525D") +
  geom_line(data = lamDF, aes(x = Year, y = Mean), colour = "#7A525D", linewidth = 0.8) +
  scale_x_continuous(limits = c(2008, 2025),
                     breaks = c(2008, 2010, 2012, 2014, 2016, 2018, 2020, 2022, 2024)) +
  scale_y_continuous(limits = c(pop_min, pop_max), breaks = c(100, 200, 300, 400, 500),
                     sec.axis = sec_axis(~ (. - s_shift) / s_fac,
                                         name = expression("Population growth rate" ~ (lambda)),
                                         breaks = c(0.6, 0.8,1.0, 1.2, 1.4))) +
  labs(x = "Year", y = "Population size") +
  theme_bw() +
  theme(legend.position = "none",
        # axis.title.y.right = element_text(color = "#7A525D"),
        # axis.text.y.right  = element_text(color = "#7A525D"),
        # axis.ticks.y.right = element_line(color = "#7A525D")
        )


# ggsave("figures/results/final!/pop&lambda_over.png", width = 18.0, height = 10.0, units = c("cm"), dpi = 600)


## Age structure ---------------------------------------------------------------

# load results
paramSamples <- readRDS('results/paramSamples.rds')

# extract proportions
pYF  <- paramSamples$t$pYF
pSA  <- paramSamples$t$pSA
pAD  <- paramSamples$t$pAD
nTOT <- paramSamples$t$nTOT

nAD   <- dim(pAD)[2]
nYear <- dim(pYF)[2]

# summarise posterior means
pYF_mean <- colMeans(pYF, na.rm = TRUE)
pSA_mean <- colMeans(pSA, na.rm = TRUE)
pAD_mean <- apply(pAD, c(2, 3), mean, na.rm = TRUE)
nTOT_mean <- colMeans(nTOT, na.rm = TRUE)

# build summary dataframe
df_YF <- tibble(Year = 1:nYear, Age = 0, Prop = pYF_mean)
df_SA <- tibble(Year = 1:nYear, Age = 1, Prop = pSA_mean)

df_AD <- expand.grid(Age = 1:nAD, Year = 1:nYear) %>%
  mutate(Prop = as.vector(pAD_mean)) %>% 
  filter(Age > 1)

df <- bind_rows(df_YF, df_SA, df_AD) %>%
  mutate(Year = Year + 2007,
         AgeGroup = case_when(Age == 0 ~ "0",
                              Age == 1 ~ "1",
                              Age %in% 2:11 ~ as.character(Age),
                              Age >= 12 ~ "12+")) %>% 
  group_by(Year, AgeGroup) %>% 
  summarise(Prop = sum(Prop), .groups = "drop") %>% 
  group_by(Year) %>%
  mutate(missing = 1 - sum(Prop),
         Prop = ifelse(AgeGroup == "12+", Prop + missing, Prop)) %>%
  ungroup() %>%
  select(-missing) %>%
  mutate(AgeGroup = factor(AgeGroup, levels = c("0", "1", as.character(2:11), "12+")),
         N = Prop * rep(nTOT_mean, each = n_distinct(AgeGroup)))

# pick colours
cols <- c("0"   = "#671313",
          "1"   = "#752915",
          setNames(colorRampPalette(c("#9E691A", "#E7D18B"))(3), as.character(2:4)),
          setNames(colorRampPalette(c("#FFF3B0", "#798C69"))(5), as.character(5:9)),
          setNames(colorRampPalette(c("#426061", "#1C3236"))(3), as.character(c(10:11, "12+"))))

# # pick colours
# cols <- c("0"   = "#CC7AA5",
#           "1"   = "#CB6723",
#           "2"   = "#AE581E",
#           setNames(colorRampPalette(c("#53A0DF", "#164C79"))(4), as.character(3:6)),
#           setNames(colorRampPalette(c("#F8C54F", "#D79909"))(3), as.character(7:9)),
#           setNames(colorRampPalette(c("#068866", "#046249"))(3), as.character(c(10:11, "12+"))))

# ribbon plot
pAGE <- df %>% 
  arrange(Year, AgeGroup) %>% # desc() to flip order
  group_by(Year) %>%
  mutate(ymin = cumsum(lag(Prop, default = 0)),
         ymax = cumsum(Prop)) %>%
  # mutate(ymin = cumsum(lag(N, default = 0)),
  # ymax = cumsum(N)) %>%
  ungroup() %>%
  ggplot(aes(x = Year, fill = AgeGroup)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), colour = NA) +
  scale_fill_manual(values = cols) +
  guides(fill = guide_legend(reverse = TRUE)) +
  scale_x_continuous(limits = c(2009, 2025),
                     breaks = c(seq(2009, 2025, by = 2))) +
  # scale_y_continuous(limits = c(0, 750)) +
  labs(x = "Year", y = "Population size", fill = "Age") +
  theme_bw(); pAGE

# ggsave("figures/results/final!/propsRibbons2.png", width = 18.0, height = 10.0, units = c("cm"), dpi = 600)

# bar plot
pAGE <- df %>%
  ggplot(aes(x = Year, y = N, fill = AgeGroup)) +
  geom_col(width = 0.8) +
  scale_fill_manual(values = cols) +
  scale_x_continuous(limits = c(2007, 2025),
                     breaks = c(2008, 2012, 2016, 2020, 2024)) +
  labs(x = "Year", y = "Proportion of the population", fill = "Age") +
  theme_bw(); pAGE

# ggsave("figures/results/final!/NsBars.png", width = 18.0, height = 10.0, units = c("cm"), dpi = 600)

# summaries to report
p0  <- paramSamples$t.mean$pYF.mean
p1  <- paramSamples$t.mean$pSA.mean
pAD <- paramSamples$t.mean$pAD.mean[,2:18]

nIter <- length(p0)

p2  <- rowSums(pAD[, 1:11, drop = FALSE])
p12 <- rowSums(pAD[, 12:17, drop = FALSE])

post <- tibble(age0 = p0, age1 = p1, age2 = p2, age12 = p12) %>% 
  mutate(total = age0 + age1 + age2 + age12,
         age12 = age12 + (1 - total)) %>%
  select(-total)

summary <- post %>%
  pivot_longer(everything(), names_to = "AgeClass", values_to = "Prop") %>%
  group_by(AgeClass) %>%
  summarise(Mean  = mean(Prop),
            Lower = quantile(Prop, 0.025),
            Upper = quantile(Prop, 0.975),
            .groups = "drop") %>%
  mutate(AgeClass = recode(AgeClass,
                           age0    = "0",
                           age1    = "1",
                           age2_11 = "2–11",
                           age12p  = "12+"))


## SIMULATED age structure -----------------------------------------------------

# VERSION 1
# Define parameters
years <- 2010:2020
ages <- c(as.character(0:11), "12+")
Nages <- length(ages)

# 1. Define baseline proportions
baseProp <- c(
  0.09, 0.07, 0.04,  # Juveniles (0, 1, 2; ~20%)
  rep(0.50 / 7, 7),  # Prime-aged (3 through 9; ~50%)
  0.10, 0.10, 0.10   # Senescent (10, 11, 12+; ~30%)
)

# 2. Population trajectory (harder crash in Year 2, smooth recovery)
totalN <- c(500, 240, 270, 310, 350, 390, 430, 460, 482, 495, 500)

# 3. Perturbation multipliers for Year 2 (heavy decimation of young and old)
pertMult <- c(
  0.10, 0.10, 0.15,  # Juveniles severely hit
  rep(0.90, 7),      # Prime-aged mostly spared
  0.15, 0.10, 0.05   # Senescent severely wiped out
)

# Build proportion matrix across the 10 years with smooth recovery
propMatrix <- matrix(0, nrow = 11, ncol = Nages)
colnames(propMatrix) <- ages

# Year 1: baseline equilibrium
propMatrix[1, ] <- baseProp

# Year 2: post-disturbance proportions
rawPert <- baseProp * pertMult
propMatrix[2, ] <- rawPert / sum(rawPert)

# Years 3 to 10: gradual linear recovery
for (t in 3:11) {
  w <- (t - 2) / 8  # progress weight from 0 to 1
  interpolated <- (1 - w) * propMatrix[2, ] + w * baseProp
  propMatrix[t, ] <- interpolated / sum(interpolated)
}

# 4. Construct the final data frame with all 13 age classes
df <- map_dfr(1:11, function(i){
  props <- propMatrix[i, ]
  totN  <- totalN[i]
  data.frame(
    Year     = years[i],
    AgeGroup = factor(ages, levels = ages),
    Prop     = props,
    N        = props * totN
  )
}) %>%
  mutate(Stage = case_when(
    AgeGroup %in% c("0", "1", "2") ~ "Juvenile",
    AgeGroup %in% as.character(3:9) ~ "Prime-aged",
    TRUE ~ "Senescent"),
    Stage = factor(Stage, levels = c("Juvenile", "Prime-aged", "Senescent")))

# 5. Group into 3 stages for plotting 
df <- df %>%
  group_by(Year, Stage) %>%
  summarise(N = sum(N), .groups = "drop")

# pick colours
# cols <- c(
#   "Juvenile"   = "#AC2020",
#   "Prime-aged" = "#E3CA78",
#   "Senescent"  = "#6F9A9B"
# )

cols <- c(
  "Juvenile"   = "#671313",
  "Prime-aged" = "#E7D18B",
  "Senescent"  = "#1C3236"
)

ggplot(df, aes(x = as.factor(Year), y = N, fill = Stage)) +
  geom_col(width = 0.8) +
  scale_fill_manual(values = cols, name = "Stage") +
  scale_x_discrete(breaks = c(2010, 2012, 2014, 2016, 2018)) +
  labs(x = "Year", y = "Population size") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    legend.position = "right"
  )

# ggsave("figures/results/final!/AgeStructEx2.png", width = 18.0, height = 10.0, units = c("cm"), dpi = 600)

# VERSION 2
# Define parameters (2010 to 2026: 17 years)
years <- 2010:2026
Nyears <- length(years)
ages <- c(as.character(0:11), "12+")
Nages <- length(ages)

# 1. Define baseline proportions
baseProp <- c(
  0.09, 0.07, 0.04,  # Juveniles (0, 1, 2; ~20%)
  rep(0.50 / 7, 7),  # Prime-aged (3 through 9; ~50%)
  0.10, 0.10, 0.10   # Senescent (10, 11, 12+; ~30%)
)

# 2. Population trajectory:
# - 2010-2020: original crash & recovery (indices 1 to 11)
# - 2021-2026: second crash & partialrecovery through 2026
totalN <- c(
  500, 240, 270, 310, 350, 390, 430, 460, 482, 495, 500, # 2010–2020
  360, 390, 415, 435, 450, 460                           # 2021–2026
)

# 3. Perturbation multipliers
# First perturbation (Year 2 / 2011): hits young & old hard
pertMult1 <- c(
  0.10, 0.10, 0.15,  
  rep(0.90, 7),      
  0.15, 0.10, 0.05   
)

# Second perturbation (Year 11 / 2020): hits prime-aged hardest
pertMult2 <- c(
  0.85, 0.85, 0.90,
  0.30, 0.35, 0.35, 0.40, 0.40, 0.45, 0.50,
  0.80, 0.75, 0.70
)

# Build proportion matrix across the 17 years
propMatrix <- matrix(0, nrow = Nyears, ncol = Nages)
colnames(propMatrix) <- ages

# Year 1 (2010): baseline equilibrium
propMatrix[1, ] <- baseProp

# Year 2 (2011): first post-disturbance proportions
rawPert1 <- baseProp * pertMult1
propMatrix[2, ] <- rawPert1 / sum(rawPert1)

# Years 3 to 11 (2012–2020): first recovery period
for (t in 3:11) {
  w <- (t - 2) / 9  
  interpolated <- (1 - w) * propMatrix[2, ] + w * baseProp
  propMatrix[t, ] <- interpolated / sum(interpolated)
}

# Year 12 (2021): second disturbance hits
rawPert2 <- propMatrix[11, ] * pertMult2
propMatrix[12, ] <- rawPert2 / sum(rawPert2)

# Years 13 to 17 (2022–2026): second partial recovery period
for (t in 13:17) {
  w <- (t - 11) / 8  # slower recovery weight
  interpolated <- (1 - w) * propMatrix[12, ] + w * baseProp
  propMatrix[t, ] <- interpolated / sum(interpolated)
}

# 4. Construct the final data frame with all age classes
df <- map_dfr(1:Nyears, function(i){
  props <- propMatrix[i, ]
  totN  <- totalN[i]
  data.frame(
    Year     = years[i],
    AgeGroup = factor(ages, levels = ages),
    Prop     = props,
    N        = props * totN
  )
}) %>%
  mutate(Stage = case_when(
    AgeGroup %in% c("0", "1", "2") ~ "Juvenile",
    AgeGroup %in% as.character(3:9) ~ "Prime-aged",
    TRUE ~ "Senescent"),
    Stage = factor(Stage, levels = c("Juvenile", "Prime-aged", "Senescent")))

# 5. Group into 3 stages for plotting  
df <- df %>%
  group_by(Year, Stage) %>%
  summarise(N = sum(N), .groups = "drop")

# pick colours
cols <- c(
  "Juvenile"   = "#963E25",
  "Prime-aged" = "#DCCC84",
  "Senescent"  = "#4D6B6C"
)

ggplot(df, aes(x = as.factor(Year), y = N, fill = Stage)) +
  geom_col(width = 0.8) +
  scale_fill_manual(values = cols, name = "Stage") +
  scale_x_discrete(breaks = c(2010, 2014, 2018, 2022, 2026)) +
  labs(x = "Year", y = "Population size") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    legend.position = "right"
  )

# ggsave("figures/results/final!/AgeStructEx3.png", width = 18.0, height = 8.0, units = c("cm"), dpi = 600)


## Dave vs Heloise's density data ----------------------------------------------

# wrangle Dave's data
source('R/wrangleData_en.R')
dave <- wrangleData_en(
  # dens.data = "data/abundanceData_Proteus.csv" # OR...
  dens.data = "data/WPNP_Methods_Results_January2026.xlsx",
  veg.data  = "data/biomass data April 2009 - July 2025_updated Feb2026.xlsx",
  wea.data  = "data/Prom_Weather_2008-2023_updated Jan2026 RB.xlsx",
  wind.data = "data/POWER_Point_Daily_20080101_20260331_10M.csv",
  obs.data  = "data/PromObs_2008-2024.xlsx",
  list.data = "data/PromlistAllNov25.xlsx")

year <- seq(2008, 2025, 1)
dave <- as.data.frame(cbind(year, dave$dens, dave$densE)) %>% 
  rename(dens = 'V2', densE = 'V3') %>% 
  mutate(data = "Dave") %>% 
  filter(year > 2008)

# wrangle Heloise's data
source('R/wrangleData_en.R')
heloise <- wrangleData_en(
  dens.data = "data/abundanceData_Proteus.csv", # OR...
  # dens.data = "data/WPNP_Methods_Results_January2026.xlsx",
  veg.data  = "data/biomass data April 2009 - July 2025_updated Feb2026.xlsx",
  wea.data  = "data/Prom_Weather_2008-2023_updated Jan2026 RB.xlsx",
  wind.data = "data/POWER_Point_Daily_20080101_20260331_10M.csv",
  obs.data  = "data/PromObs_2008-2024.xlsx",
  list.data = "data/PromlistAllNov25.xlsx")

heloise <- as.data.frame(cbind(year, heloise$dens, heloise$densE)) %>% 
  rename(dens = 'V2', densE = 'V3') %>% 
  mutate(data = "Heloise") %>% 
  filter(year > 2008)

# combine
densities <- rbind(dave, heloise)

# plot
densities %>% 
  ggplot(aes(x = year, y = dens, group = data, colour = data, fill = data)) +
  geom_ribbon(aes(ymin = dens-densE, ymax = dens+densE), alpha = 0.2, colour = NA) + 
  geom_line(linewidth = 0.8) +
  scale_x_continuous(breaks = c(2009, 2012, 2015, 2018, 2021, 2024)) +
  scale_y_continuous(limits = c(1, 10),
                     breaks = c(2, 4, 6, 8, 10)) +
  labs(x = "Year", y = "Density estimate", colour = "Data", fill = "Data") +
  theme_bw()

# ggsave("figures/densData_shrunkCIs.png", width = 18.0, height = 10.0, units = c("cm"), dpi = 600)

tmp <- as.data.frame(cbind(dave$year, dave$dens, heloise$dens)) %>% 
  rename(year = 'V1', dave = 'V2', heloise = 'V3')

library(ggpubr)
ggscatter(tmp, x = "dave", y = "heloise",
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "pearson",
          alpha = 0.2, ggtheme = theme_bw())


## Age-specific intercepts -----------------------------------------------------

postPaths <- c(
  "results/IPM_CJSen_RSen_AB_DynDens_dCJS_12_noW_25BR_S33_R22_stochV_8chains.rds",
  "results/IPM_CJSen_RSen_AB_DynDens_dCJS_12_noW_25BR_S33_R22_stochV_8chains_muFun.rds",
  "results/IPM_CJSen_RSen_AB_DynDens_dCJS_12_noW_25BR_S33_R22_stochV_8chains_muFunE.rds",
  "results/IPM_CJSen_RSen_AB_DynDens_dCJS_12_noW_25BR_S33_R22_B22_stochV_8chains_muFunE.rds"
)

modelNames <- c(
  "IPM_S33R22",
  "IPM_muFun",
  "IPM_muFunE",
  "IPM_muFunE_B22"
)

nModels <- length(modelNames)
plotFolder <- c("figures/densityChecks/muFun")

# reformat posterior samples
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

# set plotting colors
plot.cols <- paletteer_c("grDevices::Temps", nModels)

# filter posterior data for intercepts
mu.dat <- post.dat %>%
  filter(grepl("^Mu\\.[SBR]\\[\\d+\\]$", Parameter)) %>%
  mutate(
    ParamType = stringr::str_extract(Parameter, "^Mu\\.[SBR]"),
    Age = as.numeric(stringr::str_extract(Parameter, "\\d+"))
  )

# a) Survival intercepts (Mu.S)
p.Mu.S <- ggplot(subset(mu.dat, ParamType == "Mu.S")) +
  geom_violin(aes(x = factor(Age), y = Value, fill = Model, colour = Model), 
              alpha = 0.5, scale = "width", draw_quantiles = 0.5, 
              position = position_dodge(width = 0.8)) +
  ylab("Parameter value") +
  xlab("") +
  labs(title = "a) Survival intercepts") +
  scale_x_discrete(labels = c("0", 1:11, "12+")) +
  scale_fill_manual(values = plot.cols) +
  scale_colour_manual(values = plot.cols) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.title = element_text(size = 10),
        plot.margin = margin(1, 3, 1, 3)); p.Mu.S

# b) Birth rate intercepts (Mu.B)
p.Mu.B <- ggplot(subset(mu.dat, ParamType == "Mu.B")) +
  geom_violin(aes(x = factor(Age), y = Value, fill = Model, colour = Model), 
              alpha = 0.5, scale = "width", draw_quantiles = 0.5, 
              position = position_dodge(width = 0.8)) +
  ylab("Parameter value") +
  xlab("") +
  labs(title = "b) Birth rate intercepts") +
  scale_x_discrete(labels = c(2:12)) +
  scale_fill_manual(values = plot.cols) +
  scale_colour_manual(values = plot.cols) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.title = element_text(size = 10),
        plot.margin = margin(1, 3, 1, 3)); p.Mu.B

# c) Reproductive success intercepts (Mu.R)
p.Mu.R <- ggplot(subset(mu.dat, ParamType == "Mu.R")) +
  geom_violin(aes(x = factor(Age), y = Value, fill = Model, colour = Model), 
              alpha = 0.5, scale = "width", draw_quantiles = 0.5, 
              position = position_dodge(width = 0.8)) +
  ylab("Parameter value") +
  xlab("Age") +
  labs(title = "c) Survival of PYs intercepts") +
  scale_x_discrete(labels = c(2:12)) +
  scale_fill_manual(values = plot.cols) +
  scale_colour_manual(values = plot.cols) +
  theme_bw() +
  theme(legend.position = "bottom", # Keep legend only on the final bottom plot
        legend.title = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.title = element_text(size = 10),
        plot.margin = margin(1, 3, 1, 3)); p.Mu.R

# Combine with patchwork & save 
pdf(paste0(plotFolder, "/PostSummaries_MuIntercepts.pdf"), width = 10, height = 10)
print(p.Mu.S / p.Mu.B / p.Mu.R)
dev.off()

