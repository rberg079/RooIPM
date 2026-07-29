# 20 April 2026
# Plot results of IPM

## Set up ----------------------------------------------------------------------

library(tidyverse)
library(data.table)
library(patchwork)
library(scales)

nYear <- 18
nAge  <- 19

# load results
out.mcmc <- readRDS('results/IPM_CJSen_RSen_AB_DynDens_dCJS_12_noW_25BR_Dave2Covs_stochV_8chains.rds')
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

df <- data.frame(
  Year = 1:length(var),
  Mean = apply(out.mat[, var, drop = FALSE], 2, mean, na.rm = TRUE),
  Lower = apply(out.mat[, var, drop = FALSE], 2, quantile, probs = 0.025, na.rm = TRUE),
  Upper = apply(out.mat[, var, drop = FALSE], 2, quantile, probs = 0.975, na.rm = TRUE)
)

# population plot
pop <- df %>% 
  # filter(Year > 1) %>% 
  mutate(Year = Year + 2007) %>% 
  ggplot(aes(x = Year, y = Mean)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "#62414B", alpha = 0.2) + #C398B7
  geom_line(color = "#62414B", linewidth = 0.8) + #673C5B
  # scale_x_continuous(limits = c(2009, 2025),
  #                    breaks = c(2009, 2011, 2013, 2015, 2017, 2019, 2021, 2023, 2025)) +
  scale_x_continuous(limits = c(2008, 2025),
                     breaks = c(2008, 2010, 2012, 2014, 2016, 2018, 2020, 2022, 2024)) +
  scale_y_continuous(limits = c(100, 520),
                     breaks = c(100, 200, 300, 400, 500)) +
  labs(y = "Population size") +
  theme_bw() +
  theme(
        axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none"
        ); pop

# ggsave("figures/resultsDave2Covs/nTOT.jpeg", width = 18.0, height = 10.0, units = c("cm"), dpi = 600)


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
df_YF   <- data.frame(Year = 1:nYear, Group = "Young-at-foot (0 years)",
                      Mean  = apply(out.mat[, nYF_idx, drop = FALSE], 2, mean, na.rm = TRUE),
                      Lower = apply(out.mat[, nYF_idx, drop = FALSE], 2, quantile, probs = 0.025, na.rm = TRUE),
                      Upper = apply(out.mat[, nYF_idx, drop = FALSE], 2, quantile, probs = 0.975, na.rm = TRUE))

df_SA   <- data.frame(Year = 1:nYear, Group = "Subadults (1 year)",
                      Mean  = apply(out.mat[, nSA_idx, drop = FALSE], 2, mean, na.rm = TRUE),
                      Lower = apply(out.mat[, nSA_idx, drop = FALSE], 2, quantile, probs = 0.025, na.rm = TRUE),
                      Upper = apply(out.mat[, nSA_idx, drop = FALSE], 2, quantile, probs = 0.975, na.rm = TRUE))

df_2to9 <- data.frame(Year = 1:nYear, Group = "Adults (2-9 years)",
                      Mean  = apply(mat_2to9, 2, mean, na.rm = TRUE),
                      Lower = apply(mat_2to9, 2, quantile, probs = 0.025, na.rm = TRUE),
                      Upper = apply(mat_2to9, 2, quantile, probs = 0.975, na.rm = TRUE))

df_10up <- data.frame(Year = 1:nYear, Group = "Adults (10+ years)",
                      Mean  = apply(mat_10up, 2, mean, na.rm = TRUE),
                      Lower = apply(mat_10up, 2, quantile, probs = 0.025, na.rm = TRUE),
                      Upper = apply(mat_10up, 2, quantile, probs = 0.975, na.rm = TRUE))

df_TOT  <- data.frame(Year = 1:nYear, Group = "Total female population",
                      Mean  = apply(out.mat[, nTOT_idx, drop = FALSE], 2, mean, na.rm = TRUE),
                      Lower = apply(out.mat[, nTOT_idx, drop = FALSE], 2, quantile, probs = 0.025, na.rm = TRUE),
                      Upper = apply(out.mat[, nTOT_idx, drop = FALSE], 2, quantile, probs = 0.975, na.rm = TRUE))

# combine everything
df <- rbind(df_YF, df_SA, df_2to9, df_10up, df_TOT)

# pick colours
# cols <- c(
#   "YF"       = "#F9C74F", # Yellow 
#   "SA"       = "#F3722C", # Orange
#   "AD (2-9)" = "#D62828", # Red
#   "AD (10+)" = "#8E5EA2", # Purple (borrowed from survival script)
#   "TOT"      = "#62414B"  # Dark purple/brown
# )

cols <- c(
  "Young-at-foot (0 years)" = "#8E5EA2",
  "Subadults (1 year)"      = "#277DA1",
  "Adults (2-9 years)"      = "#F9C74F",
  "Adults (10+ years)"      = "#D62828",
  "Total female population" = "#62414B"
)

# population plot
popA <- df %>% 
  # filter(Year > 1) %>% 
  mutate(Year = Year + 2007,
         Group = factor(Group,
                        levels = c("Young-at-foot (0 years)",
                                   "Subadults (1 year)",
                                   "Adults (2-9 years)",
                                   "Adults (10+ years)",
                                   "Total female population"))) %>% 
  ggplot(aes(x = Year, y = Mean, group = Group, colour = Group, fill = Group)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.2, colour = NA) + 
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = cols) +
  scale_fill_manual(values = cols) +
  scale_x_continuous(limits = c(2008, 2025),
                     breaks = c(2008, 2010, 2012, 2014, 2016, 2018, 2020, 2022, 2024)) +
  scale_y_continuous(limits = c(0, 520), breaks = scales::pretty_breaks()) +
  labs(y = "Population size", colour = "Age group", fill = "Age group") +
  theme_bw(); popA

# ggsave("figures/resultsDave2Covs/allNs.jpeg", width = 18.0, height = 10.0, units = c("cm"), dpi = 600)


## Survival --------------------------------------------------------------------

# indices
S_idx  <- grep("^S\\[", colnames(out.mat))

# build summary dataframe
df <- expand.grid(Age = 1:13, Year = 1:17)
df$Mean <- apply(out.mat[, S_idx, drop = FALSE], 2, mean, na.rm = TRUE)
df$Lower <- apply(out.mat[, S_idx, drop = FALSE], 2, quantile, probs = 0.025, na.rm = TRUE)
df$Upper <- apply(out.mat[, S_idx, drop = FALSE], 2, quantile, probs = 0.975, na.rm = TRUE)

# pick colours
# cols <- c(
#   "0"  = "#CC6673",
#   "1"  = "#AB2136",
#   "2"  = "#DE612B",
#   "4"  = "#F4A261",
#   "6"  = "#E9C46A",
#   "8"  = "#2A9D8F",
#   "10" = "#264653",
#   "12+" = "#936271"
# )

# cols <- c(
#   "0"  = "#8E5EA2",
#   "1"  = "#277DA1",
#   "2"  = "#00A896",
#   "4"  = "#F9C74F",
#   "6"  = "#F8A23A",
#   "8"  = "#F3722C",
#   "10" = "#D62828",
#   "12+" = "#D96C9D"
# )

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

# ggsave("figures/resultsDave2Covs/survival.jpeg", width = 18.0, height = 10.0, units = c("cm"), dpi = 600)


## Reproductive success --------------------------------------------------------

# indices
# Bt_idx  <- grep("^Bt\\[", colnames(out.mat))
BR_idx  <- grep("^BR\\[", colnames(out.mat))
sPY_idx <- grep("^sPY\\[", colnames(out.mat))

# extract matrices
# Bt  <- out.mat[, Bt_idx,  drop = FALSE]
# Bt  <- Bt[, rep(1:ncol(Bt), each = 19)]
BR  <- out.mat[, BR_idx, drop = FALSE]
sPY <- out.mat[, sPY_idx, drop = FALSE]

# compute reproductive output
# R <- 0.5 * Bt * sPY
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

# # summarise by age & year
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

# pick colours
# cols <- c(
#   "2"  = "#DE612B",
#   "4"  = "#F4A261",
#   "6"  = "#E9C46A",
#   "8"  = "#2A9D8F",
#   "10" = "#264653",
#   "12+" = "#936271"
# )

# cols <- c(
#   "2"  = "#00A896",
#   "4"  = "#F9C74F",
#   "6"  = "#F8A23A",
#   "8"  = "#F3722C",
#   "10" = "#D62828",
#   "12+" = "#D96C9D"
# )

cols <- c(
  "2"   = "#009E73",
  "4"   = "#F0E442",
  "6"   = "#E69F00",
  "8"   = "#D55E00",
  "10"  = "#CC3311",
  "12+" = "#CC79A7"
)


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
  geom_line(linewidth = 0.8, show.legend = F) +
  scale_colour_manual(values = cols) +
  scale_fill_manual(values = cols) +
  scale_linetype_manual(values = lts) +
  scale_x_continuous(limits = c(2008, 2024),
                     breaks = c(2008, 2010, 2012, 2014, 2016, 2018, 2020, 2022, 2024)) +
  # scale_y_continuous(limits = c(0, 1), breaks = pretty_breaks()) +
  scale_y_continuous(limits = c(0, 0.5),
                     breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5)) +
  labs(x = "Year", y = "Reproductive success", colour = "Age", fill = "Age") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_blank()
  ); rs

# ggsave("figures/resultsDave2Covs/PYsurv.jpeg", width = 18.0, height = 10.0, units = c("cm"), dpi = 600)

# combine with survival plot
# surv / rs

# combine birth rate & PY survival
(br / rs) +
  plot_layout(guides = "collect") +
  theme(legend.position = "right")

# ggsave("figures/resultsDave2Covs/birth&PYsurv2.jpeg", width = 18.0, height = 18.0, units = c("cm"), dpi = 600)

# ...& population size
(surv / rs / pop) +
  plot_layout(guides = "collect") +
  theme(legend.position = "right")

# ggsave("figures/resultsDave2Covs/surv&rs&pop.jpeg", width = 18.0, height = 22.0, units = c("cm"), dpi = 600)


## Covariate effects -----------------------------------------------------------

inv_logit <- function(x) 1 / (1 + exp(-x))

keep_ages <- c(1, 11, 12, 13)
age_labels <- c("0", "10", "11", "12+")

# extract betas
bDy <- out.mat[, grep("^BetaD\\.Sy", colnames(out.mat))]
bDo <- out.mat[, grep("^BetaD\\.So", colnames(out.mat))]
bVy <- out.mat[, grep("^BetaV\\.Sy", colnames(out.mat))]
bVo <- out.mat[, grep("^BetaV\\.So", colnames(out.mat))]

bDr <- out.mat[, grep("^BetaD\\.R", colnames(out.mat))]

# extract baseline intercepts
muS  <- paste0("Mu.S[", keep_ages, "]")
base <- out.mat[, muS, drop = FALSE]

muR <- grep("^Mu\\.R\\[", colnames(out.mat))
baseR <- rowMeans(out.mat[, muR, drop = FALSE])

# covariate sequence
x <- seq(-2, 2, length.out = 50)

# build summary dfs
make_df <- function(betas_young, betas_old, covariate) {
  results <- list()
  
  for (i in seq_along(keep_ages)) {
    age <- keep_ages[i]
    label <- age_labels[i]
    intercepts <- base[, i]
    
    # assign the right beta
    if (age == 1) {
      betas <- betas_young
    } else {
      betas <- betas_old
    }
    
    for (xx in x) {
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
make_py_df <- function(betas, covariate) {
  results <- list()
  for (xx in x) {
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
  make_df(bDy, bDo, "Population density"),
  make_df(bVy, bVo, "Forage availability"),
  make_py_df(bDr, "Population density")
)

# pick colours
# cols <- c(
#   "0"   = "#CC6673",
#   "10"  = "#264653",
#   "11"  = "#5D5462",
#   "12+" = "#936271"
# )
# 
# cols <- c(
#   "0"   = "#8E5EA2",
#   "10"  = "#D62828",
#   "11"  = "#D84A63",
#   "12+" = "#D96C9D"
# )
# 
# cols <- c(
#   "0"   = "#714A82",
#   "10"  = "#9B1C1C",
#   "11"  = "#D43552",
#   "12+" = "#D96D9E"
# )

# cols <- c(
#   "PY"  = "#BFA4CB",
#   "0"   = "#7C528E",
#   "10"  = "#BE2323",
#   "11"  = "#696969",
#   "12+" = "#D96C9D"
# )

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
  # scale_y_continuous(limits = c(0, 1), breaks = c(0.2, 0.4, 0.6, 0.8, 1.0)) +
  labs(x = "Scaled covariate value", y = "Survival",
       colour = "Age", fill = "Age") +
  guides(linetype = "none") +
  theme_bw()

# ggsave("figures/resultsDave2Covs/coveffects2.jpeg", width = 18.0, height = 10.0, units = c("cm"), dpi = 600)


## Covariate values ------------------------------------------------------------

# indices
D_idx  <- grep("^D_dens\\.true", colnames(out.mat))[1:17]
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
  Year = 1:17, 
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
  "Population density"  = "#47404F",
  "Forage availability" = "#7D9570"
)

fill_cols <- c(
  "Population density"  = adjustcolor("#47404F", alpha.f = 0.15),
  "Forage availability" = adjustcolor("#7D9570", alpha.f = 0.20)
)

# plot
covs <- df %>%
  filter(Year > 1) %>% 
  mutate(Year = Year + 2007) %>%
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

# ggsave("figures/resultsDave2Covs/covsVStime2.jpeg", width = 18.0, height = 10.0, units = c("cm"), dpi = 600)

# combine with survival & reproductive output
(surv / rs / covs) +
  plot_layout(guides = "collect") +
  theme(legend.position = "right")

# ggsave("figures/resultsDave2Covs/surv&rs&covs4.jpeg", width = 20.0, height = 22.0, units = c("cm"), dpi = 600)

# plot observed vs true density
myData$D_densE[1] <- NA

D_idx_all <- grep("^D_dens\\.true", colnames(out.mat))[1:18]

true_dens_mean <- apply(out.mat[, D_idx_all, drop = FALSE], 2, mean, na.rm = TRUE)
true_dens_lci  <- apply(out.mat[, D_idx_all, drop = FALSE], 2, quantile, probs = 0.025, na.rm = TRUE)
true_dens_uci  <- apply(out.mat[, D_idx_all, drop = FALSE], 2, quantile, probs = 0.975, na.rm = TRUE)

# combine & calculate 95% CrIs for observed density 
# (D_densE is standard error: Mean +/- 1.96 * SE)
df_obs <- data.frame(
  Year = (1:18) + 2007,
  True_Mean = true_dens_mean,
  True_LCI  = true_dens_lci,
  True_UCI  = true_dens_uci,
  Obs_Mean  = myData$D_dens,
  Obs_SE    = myData$D_densE
) %>%
  mutate(
    Obs_LCI = Obs_Mean - (1.96 * Obs_SE),
    Obs_UCI = Obs_Mean + (1.96 * Obs_SE)
  )

# plot
obs_vs_true <- df_obs %>%
  filter(Year > 2008) %>%
  ggplot(aes(x = Year)) +
  
  geom_ribbon(aes(ymin = Obs_LCI, ymax = Obs_UCI, fill = "Observed"), alpha = 0.15) +
  geom_line(aes(y = Obs_Mean, colour = "Observed"), linewidth = 0.8, alpha = 0.8) +
  # geom_errorbar(aes(ymin = Obs_LCI, ymax = Obs_UCI, colour = "Observed"), width = 0.3) +
  # geom_point(aes(y = Obs_Mean, colour = "Observed"), size = 2) +
  
  geom_ribbon(aes(ymin = True_LCI, ymax = True_UCI, fill = "Estimated (IPM)"), alpha = 0.2) +
  geom_line(aes(y = True_Mean, colour = "Estimated (IPM)"), linewidth = 0.8) +
  scale_colour_manual(values = c("Estimated (IPM)" = "#335B5B", "Observed" = "#D55E00")) +
  scale_fill_manual(values = c("Estimated (IPM)" = "#335B5B", "Observed" = "#D55E00")) +
  scale_x_continuous(breaks = seq(2008, 2024, by = 2)) +
  labs(x = "Year", y = "Population density", colour = "Source", fill = "Source") +
  theme_bw(); obs_vs_true

# ggsave("figures/resultsDave2Covs/obsVStrueDens.jpeg", width = 18.0, height = 10.0, units = c("cm"), dpi = 600)


## Lambda ----------------------------------------------------------------------

# load results
paramSamples <- readRDS('results/paramSamples.rds')

# extract lambda
lambda_mat <- paramSamples$t$lambda

# summary stats
mean <- apply(lambda_mat, 1, mean)
quantile(mean, c(0.025, 0.5, 0.975))

# build summary dataframe
df <- data.frame(Year  = 1:17,
                 Mean  = apply(lambda_mat, 2, mean),
                 Lower = apply(lambda_mat, 2, quantile, 0.025),
                 Upper = apply(lambda_mat, 2, quantile, 0.975))

lambda <- df %>%
  mutate(Year = Year + 2007) %>%
  ggplot(aes(x = Year, y = Mean)) +
  geom_hline(yintercept = 1, colour = "grey60") +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.2, fill = "#8D6B48") +
  geom_line(colour = "#8D6B48", linewidth = 0.8) +
  scale_x_continuous(limits = c(2008, 2025),
                     breaks = c(2008, 2010, 2012, 2014, 2016, 2018, 2020, 2022, 2024)) +
  # scale_x_continuous(limits = c(2008, 2025),
  #                    breaks = c(2009, 2011, 2013, 2015, 2017, 2019, 2021, 2023, 2025)) +
  scale_y_continuous(limits = c(NA, 1.25),
                     breaks = c(0.6, 0.8, 1.0, 1.2)) +
  labs(x = "Year", y = expression("Population growth rate" ~ (lambda))) +
  theme_bw(); lambda

# ggsave("figures/resultsDave2Covs/lambda.jpeg", width = 18.0, height = 10.0, units = c("cm"), dpi = 600)

# combine with pop size
pop / lambda
pAGE / lambda

# ggsave("figures/resultsDave2Covs/pop&lambda.jpeg", width = 18.0, height = 18.0, units = c("cm"), dpi = 600)


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

# colour palette
# cols <- c("0"   = "#E8D6CB",
#           "1"   = "#D0ADA7",
#           setNames(colorRampPalette(c("#AD6A6C", "#5D2E46"))(10), as.character(2:11)),
#           "12+" = "#371B29")
# 
# cols <- c("0"   = "#C9CBA3",
#           "1"   = "#FFE1A8",
#           setNames(colorRampPalette(c("#EDA297", "#DD5540"))(8), as.character(2:9)),
#           setNames(colorRampPalette(c("#854752", "#5D3239"))(3), as.character(c("10", "11", "12+"))))
# 
# cols <- c("0"   = "#E3AD78",
#           "1"   = "#EFD6AC",
#           setNames(colorRampPalette(c("#A0A290", "#4D4F40"))(8), as.character(2:9)),
#           setNames(colorRampPalette(c("#695C70", "#4A404F"))(3), as.character(c("10", "11", "12+"))))
# 
# cols <- c("0"   = "#CA562C",
#           "1"   = "#DE8A5A",
#           setNames(colorRampPalette(c("#FAF5DB", "#F2D4A4"))(3), as.character(2:4)),
#           setNames(colorRampPalette(c("#B5B991", "#657359"))(5), as.character(5:9)),
#           setNames(colorRampPalette(c("#324935", "#19241A"))(3), as.character(c(10:11, "12+"))))

cols <- c("0"   = "#671313",
          "1"   = "#752915",
          setNames(colorRampPalette(c("#9E691A", "#E7D18B"))(3), as.character(2:4)),
          setNames(colorRampPalette(c("#FFF3B0", "#798C69"))(5), as.character(5:9)),
          setNames(colorRampPalette(c("#426061", "#1C3236"))(3), as.character(c(10:11, "12+"))))

# cols <- paletteer::paletteer_c("grDevices::Temps", 7)

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

# ggsave("figures/resultsDave2Covs/propsRibbons.jpeg", width = 18.0, height = 10.0, units = c("cm"), dpi = 600)

# bar plot
pAGE <- df %>%
  ggplot(aes(x = Year, y = N, fill = AgeGroup)) +
  geom_col(width = 0.8) +
  scale_fill_manual(values = cols) +
  scale_x_continuous(limits = c(2009, 2024),
                     breaks = c(2010, 2012, 2014, 2016, 2018, 2020, 2022, 2024)) +
  labs(x = "Year", y = "Proportion of the population", fill = "Age") +
  theme_bw(); pAGE

# ggsave("figures/resultsDave2Covs/NsBars.jpeg", width = 18.0, height = 10.0, units = c("cm"), dpi = 600)

# summaries of t.mean to report
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


## Sensitivity Analysis --------------------------------------------------------

library(tidyverse)
library(data.table)
library(patchwork)
library(scales)

oneProp <- TRUE # whether age structure should be summed

sensitivities <- readRDS('results/sensitivities.rds')
sens_samples <- sensitivities$sensitivity$samples 
nSamples <- length(sens_samples$sens.sYF)
nAge <- 19

# function to pull age-specific columns & calculate sums
processMatrix <- function(mat, prefix) {
  df <- as.data.frame(mat)
  
  # name age-specific columns
  colnames(df) <- paste0(prefix, "_", 1:ncol(df))
  
  # calculate & append grouped sums
  df[[paste0(prefix, "_all")]]  <- rowSums(mat)
  df[[paste0(prefix, "_2to9")]] <- rowSums(mat[, 2:9, drop = FALSE])
  df[[paste0(prefix, "_10up")]] <- rowSums(mat[, 10:ncol(mat), drop = FALSE])
  
  return(df)
}

# 1D variables
df_1D <- data.frame(
  sYF = sens_samples$sens.sYF,
  sSA = sens_samples$sens.sSA,
  pYF = sens_samples$sens.pYF,
  pSA = sens_samples$sens.pSA
)

# 2D variables
df_BR  <- processMatrix(sens_samples$sens.BR, "BR")
df_sPY <- processMatrix(sens_samples$sens.sPY, "sPY")
df_sAD <- processMatrix(sens_samples$sens.sAD, "sAD")
df_pAD <- processMatrix(sens_samples$sens.pAD, "pAD")

# bind everything together
wideData <- bind_cols(df_1D, df_BR, df_sPY, df_sAD, df_pAD)
wideData$draw <- 1:nrow(wideData)

# pivot to long format
sensData <- wideData %>% 
  pivot_longer(-draw, names_to = "Variable", values_to = "Sensitivity") %>% 
  mutate(Parameter = str_split_fixed(Variable, "_", 2)[, 1],
         Age = str_split_fixed(Variable, "_", 2)[, 2]) %>% 
  mutate(Age = ifelse(Age == "", NA, Age))

# format for plotting
if(oneProp) {
  tmp <- sensData %>% 
    filter(Variable %in% c("pYF", "pSA", "pAD_2to9", "pAD_10up")) %>% 
    group_by(draw) %>% 
    summarise(Sensitivity = sum(Sensitivity), .groups = "drop") %>% 
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
  
  remove(tmp)
  
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

# pick colours
plotData$type <- factor(plotData$type, levels = types)
plot.colours <- paletteer::paletteer_c("grDevices::Temps", length(unique(plotData$type)))

p.sum <- ggplot(plotData, aes(x = type, y = Sensitivity, group = type)) +
  geom_violin(aes(fill = type, colour = type), alpha = 0.5, scale = "width", draw_quantiles = 0.5) +
  geom_hline(yintercept = 0, colour = "grey70", linetype = "dashed") +
  ylab("Sensitivity") +
  xlab("") +
  scale_fill_manual(values = plot.colours) +
  scale_colour_manual(values = plot.colours) +
  scale_x_discrete(labels = names) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.title = element_text(size = 10)); p.sum

# ggsave("figures/resultsDave2Covs/SENSsum.jpeg", width = 20.0, height = 12.0, units = c("cm"), dpi = 600)

# survival panel
S.colours <- c(plot.colours[3:4], rep(plot.colours[5], 8), rep(plot.colours[6], 10))
names(S.colours) <- c("sYF", "sSA", paste0("sAD_", 2:nAge))

p.S <- ggplot(subset(sensData, Variable %in% c("sYF", "sSA", paste0("sAD_", 2:nAge)))) +
  geom_violin(aes(x = factor(Variable, levels = c("sYF", "sSA", paste0("sAD_", 2:nAge))),
                  y = Sensitivity, fill = Variable), alpha = 0.5, scale = "width", draw_quantiles = 0.5) +
  geom_hline(yintercept = 0, colour = "grey70", linetype = "dashed") +
  ylab("Sensitivity") +
  xlab("") +
  labs(title = "a) Survival") +
  scale_x_discrete(labels = expression(
    italic(S)[0], italic(S)[1],
    italic(S)[2], italic(S)[3], italic(S)[4], italic(S)[5], italic(S)[6],
    italic(S)[7], italic(S)[8], italic(S)[9], italic(S)[10], italic(S)[11],
    italic(S)[12], italic(S)[13], italic(S)[14], italic(S)[15], italic(S)[16],
    italic(S)[17], italic(S)[18], italic(S)[19]
  )) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  scale_fill_manual(values = S.colours) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.title = element_text(size = 10),
        plot.margin = margin(1, 3, 1, 3)); p.S

# birth rate panel
B.colours <- c(rep(plot.colours[1], 18))
names(B.colours) <- c(paste0("BR_", 2:nAge))

p.B <- ggplot(subset(sensData, Variable %in% c(paste0("BR_", 2:nAge)))) +
  geom_violin(aes(x = factor(Variable, levels = c(paste0("BR_", 2:nAge))),
                  y = Sensitivity, fill = Variable), alpha = 0.5, scale = "width", draw_quantiles = 0.5) +
  geom_hline(yintercept = 0, colour = "grey70", linetype = "dashed") +
  ylab("Sensitivity") +
  xlab("") +
  labs(title = "b) Birth rate") +
  scale_x_discrete(labels = expression(
    italic(B)[2], italic(B)[3], italic(B)[4], italic(B)[5], italic(B)[6],
    italic(B)[7], italic(B)[8], italic(B)[9], italic(B)[10], italic(B)[11],
    italic(B)[12], italic(B)[13], italic(B)[14], italic(B)[15], italic(B)[16],
    italic(B)[17], italic(B)[18], italic(B)[19]
  )) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  scale_fill_manual(values = B.colours) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.title = element_text(size = 10),
        plot.margin = margin(1, 3, 1, 3)); p.B

# reproductive success panel
R.colours <- c(rep(plot.colours[2], 18))
names(R.colours) <- c(paste0("sPY_", 2:nAge))

p.R <- ggplot(subset(sensData, Variable %in% c(paste0("sPY_", 2:nAge)))) +
  geom_violin(aes(x = factor(Variable, levels = c(paste0("sPY_", 2:nAge))),
                  y = Sensitivity, fill = Variable), alpha = 0.5, scale = "width", draw_quantiles = 0.5) +
  geom_hline(yintercept = 0, colour = "grey70", linetype = "dashed") +
  ylab("Sensitivity") +
  xlab("") +
  labs(title = "c) Survival of pouch young") +
  scale_x_discrete(labels = expression(
    italic(SP)[2], italic(SP)[3], italic(SP)[4], italic(SP)[5], italic(SP)[6],
    italic(SP)[7], italic(SP)[8], italic(SP)[9], italic(SP)[10], italic(SP)[11],
    italic(SP)[12], italic(SP)[13], italic(SP)[14], italic(SP)[15], italic(SP)[16],
    italic(SP)[17], italic(SP)[18], italic(SP)[19]
  )) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  scale_fill_manual(values = R.colours) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.title = element_text(size = 10),
        plot.margin = margin(1, 3, 1, 3)); p.R

# population structure panel
P.colours <- c(rep(plot.colours[7], 20))
names(P.colours) <- c("pYF", "pSA", paste0("pAD_", 2:nAge))

p.P <- ggplot(subset(sensData, Variable %in% c("pYF", "pSA", paste0("pAD_", 2:nAge)))) +
  geom_violin(aes(x = factor(Variable, levels = c("pYF", "pSA", paste0("pAD_", 2:nAge))),
                  y = Sensitivity, fill = Variable), alpha = 0.5, scale = "width", draw_quantiles = 0.5) +
  geom_hline(yintercept = 0, colour = "grey70", linetype = "dashed") +
  ylab("Sensitivity") +
  xlab("") +
  labs(title = "d) Population proportions") +
  scale_x_discrete(labels = expression(
    italic(P)[0], italic(P)[1],
    italic(P)[2], italic(P)[3], italic(P)[4], italic(P)[5], italic(P)[6],
    italic(P)[7], italic(P)[8], italic(P)[9], italic(P)[10], italic(P)[11],
    italic(P)[12], italic(P)[13], italic(P)[14], italic(P)[15], italic(P)[16],
    italic(P)[17], italic(P)[18], italic(P)[19]
  )) +
  scale_y_continuous(limits = c(0, 1.5), expand = expansion(mult = c(0, 0.02))) +
  scale_fill_manual(values = P.colours) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.title = element_text(size = 10),
        plot.margin = margin(1, 3, 1, 3)); p.P

# combine panels
# ((p.S + labs(tag = "a)")) / (p.B + labs(tag = "b)")) / (p.R + labs(tag = "c)")) / (p.P + labs(tag = "d)")))
p.S / p.B / p.R / p.P

# ggsave("figures/resultsDave2Covs/SENSage2.jpeg", width = 20.0, height = 24.0, units = c("cm"), dpi = 600)

# summaries to report
sensSummary <- sensData %>%
  group_by(Variable) %>%
  summarise(Mean  = mean(Sensitivity, na.rm = TRUE),
            Lower = quantile(Sensitivity, 0.025, na.rm = TRUE),
            Upper = quantile(Sensitivity, 0.975, na.rm = TRUE),
            .groups = "drop")


## Elasticity Analysis ---------------------------------------------------------

library(tidyverse)
library(data.table)
library(patchwork)
library(scales)

oneProp <- TRUE # whether age structure should be summed

sensitivities <- readRDS('results/sensitivities.rds')
elas_samples <- sensitivities$elasticity$samples 
nSamples <- length(elas_samples$elas.sYF)
nAge <- 19

# function to pull age-specific columns & calculate sums
processMatrix <- function(mat, prefix) {
  df <- as.data.frame(mat)
  
  # name age-specific columns
  colnames(df) <- paste0(prefix, "_", 1:ncol(df))
  
  # calculate & append grouped sums
  df[[paste0(prefix, "_all")]]  <- rowSums(mat)
  df[[paste0(prefix, "_2to9")]] <- rowSums(mat[, 2:9, drop = FALSE])
  df[[paste0(prefix, "_10up")]] <- rowSums(mat[, 10:ncol(mat), drop = FALSE])
  
  return(df)
}

# 1D variables
df_1D <- data.frame(
  sYF = elas_samples$elas.sYF,
  sSA = elas_samples$elas.sSA,
  pYF = elas_samples$elas.pYF,
  pSA = elas_samples$elas.pSA
)

# 2D variables
df_BR  <- processMatrix(elas_samples$elas.BR, "BR")
df_sPY <- processMatrix(elas_samples$elas.sPY, "sPY")
df_sAD <- processMatrix(elas_samples$elas.sAD, "sAD")
df_pAD <- processMatrix(elas_samples$elas.pAD, "pAD")

# bind everything together
wideData <- bind_cols(df_1D, df_BR, df_sPY, df_sAD, df_pAD)
wideData$draw <- 1:nrow(wideData)

# pivot to long format
elasData <- wideData %>% 
  pivot_longer(-draw, names_to = "Variable", values_to = "Elasticity") %>% 
  mutate(Parameter = str_split_fixed(Variable, "_", 2)[, 1],
         Age = str_split_fixed(Variable, "_", 2)[, 2]) %>% 
  mutate(Age = ifelse(Age == "", NA, Age))

# format for plotting
if(oneProp) {
  tmp <- elasData %>% 
    filter(Variable %in% c("pYF", "pSA", "pAD_2to9", "pAD_10up")) %>% 
    group_by(draw) %>% 
    summarise(Elasticity = sum(Elasticity), .groups = "drop") %>% 
    mutate(Variable = "p_all", type = "Population structure")
  
  plotData <- elasData %>% 
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
  
  remove(tmp)
  
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
  plotData <- elasData %>% 
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

# pick colours
plotData$type <- factor(plotData$type, levels = types)
plot.colours <- paletteer::paletteer_c("grDevices::Temps", length(unique(plotData$type)))

e.sum <- ggplot(plotData, aes(x = type, y = Elasticity, group = type)) +
  geom_violin(aes(fill = type, colour = type), alpha = 0.5, scale = "width", draw_quantiles = 0.5) +
  geom_hline(yintercept = 0, colour = "grey70", linetype = "dashed") +
  ylab("Elasticity") +
  xlab("") +
  scale_fill_manual(values = plot.colours) +
  scale_colour_manual(values = plot.colours) +
  scale_x_discrete(labels = names) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.title = element_text(size = 10)); e.sum

# ggsave("figures/resultsDave2Covs/ELASsum.jpeg", width = 20.0, height = 12.0, units = c("cm"), dpi = 600)

# survival panel
S.colours <- c(plot.colours[3:4], rep(plot.colours[5], 8), rep(plot.colours[6], 10))
names(S.colours) <- c("sYF", "sSA", paste0("sAD_", 2:nAge))

p.S <- ggplot(subset(elasData, Variable %in% c("sYF", "sSA", paste0("sAD_", 2:nAge)))) +
  geom_violin(aes(x = factor(Variable, levels = c("sYF", "sSA", paste0("sAD_", 2:nAge))),
                  y = Elasticity, fill = Variable), alpha = 0.5, scale = "width", draw_quantiles = 0.5) +
  geom_hline(yintercept = 0, colour = "grey70", linetype = "dashed") +
  ylab("Elasticity") +
  xlab("") +
  labs(title = "a) Survival") +
  scale_x_discrete(labels = expression(
    italic(S)[0], italic(S)[1],
    italic(S)[2], italic(S)[3], italic(S)[4], italic(S)[5], italic(S)[6],
    italic(S)[7], italic(S)[8], italic(S)[9], italic(S)[10], italic(S)[11],
    italic(S)[12], italic(S)[13], italic(S)[14], italic(S)[15], italic(S)[16],
    italic(S)[17], italic(S)[18], italic(S)[19]
  )) +
  scale_y_continuous(limits = c(0, 0.15), expand = expansion(mult = c(0, 0.02))) +
  scale_fill_manual(values = S.colours) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.title = element_text(size = 10),
        plot.margin = margin(1, 3, 1, 3)); p.S

# birth rate panel
B.colours <- c(rep(plot.colours[1], 18))
names(B.colours) <- c(paste0("BR_", 2:nAge))

p.B <- ggplot(subset(elasData, Variable %in% c(paste0("BR_", 2:nAge)))) +
  geom_violin(aes(x = factor(Variable, levels = c(paste0("BR_", 2:nAge))),
                  y = Elasticity, fill = Variable), alpha = 0.5, scale = "width", draw_quantiles = 0.5) +
  geom_hline(yintercept = 0, colour = "grey70", linetype = "dashed") +
  ylab("Elasticity") +
  xlab("") +
  labs(title = "b) Birth rate") +
  scale_x_discrete(labels = expression(
    italic(B)[2], italic(B)[3], italic(B)[4], italic(B)[5], italic(B)[6],
    italic(B)[7], italic(B)[8], italic(B)[9], italic(B)[10], italic(B)[11],
    italic(B)[12], italic(B)[13], italic(B)[14], italic(B)[15], italic(B)[16],
    italic(B)[17], italic(B)[18], italic(B)[19]
  )) +
  scale_y_continuous(limits = c(0, 0.03), expand = expansion(mult = c(0, 0.02))) +
  scale_fill_manual(values = B.colours) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.title = element_text(size = 10),
        plot.margin = margin(1, 3, 1, 3)); p.B

# reproductive success panel
R.colours <- c(rep(plot.colours[2], 18))
names(R.colours) <- c(paste0("sPY_", 2:nAge))

p.R <- ggplot(subset(elasData, Variable %in% c(paste0("sPY_", 2:nAge)))) +
  geom_violin(aes(x = factor(Variable, levels = c(paste0("sPY_", 2:nAge))),
                  y = Elasticity, fill = Variable), alpha = 0.5, scale = "width", draw_quantiles = 0.5) +
  geom_hline(yintercept = 0, colour = "grey70", linetype = "dashed") +
  ylab("Elasticity") +
  xlab("") +
  labs(title = "c) Survival of pouch young") +
  scale_x_discrete(labels = expression(
    italic(SP)[2], italic(SP)[3], italic(SP)[4], italic(SP)[5], italic(SP)[6],
    italic(SP)[7], italic(SP)[8], italic(SP)[9], italic(SP)[10], italic(SP)[11],
    italic(SP)[12], italic(SP)[13], italic(SP)[14], italic(SP)[15], italic(SP)[16],
    italic(SP)[17], italic(SP)[18], italic(SP)[19]
  )) +
  scale_y_continuous(limits = c(0, 0.03), expand = expansion(mult = c(0, 0.02))) +
  scale_fill_manual(values = R.colours) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.title = element_text(size = 10),
        plot.margin = margin(1, 3, 1, 3)); p.R

# population structure panel
P.colours <- c(rep(plot.colours[7], 20))
names(P.colours) <- c("pYF", "pSA", paste0("pAD_", 2:nAge))

p.P <- ggplot(subset(elasData, Variable %in% c("pYF", "pSA", paste0("pAD_", 2:nAge)))) +
  geom_violin(aes(x = factor(Variable, levels = c("pYF", "pSA", paste0("pAD_", 2:nAge))),
                  y = Elasticity, fill = Variable), alpha = 0.5, scale = "width", draw_quantiles = 0.5) +
  geom_hline(yintercept = 0, colour = "grey70", linetype = "dashed") +
  ylab("Elasticity") +
  xlab("") +
  labs(title = "d) Population proportions") +
  scale_x_discrete(labels = expression(
    italic(P)[0], italic(P)[1],
    italic(P)[2], italic(P)[3], italic(P)[4], italic(P)[5], italic(P)[6],
    italic(P)[7], italic(P)[8], italic(P)[9], italic(P)[10], italic(P)[11],
    italic(P)[12], italic(P)[13], italic(P)[14], italic(P)[15], italic(P)[16],
    italic(P)[17], italic(P)[18], italic(P)[19]
  )) +
  scale_y_continuous(limits = c(0, 0.15),expand = expansion(mult = c(0, 0.02))) +
  scale_fill_manual(values = P.colours) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.title = element_text(size = 10),
        plot.margin = margin(1, 3, 1, 3)); p.P

# combine panels
# ((p.S + labs(tag = "a)")) / (p.B + labs(tag = "b)")) / (p.R + labs(tag = "c)")) / (p.P + labs(tag = "d)")))
p.S / p.B / p.R / p.P

# ggsave("figures/resultsDave2Covs/ELASage2.jpeg", width = 20.0, height = 24.0, units = c("cm"), dpi = 600)

# summaries to report
elasSummary <- plotData %>%
  group_by(Variable) %>%
  summarise(Mean  = mean(Elasticity, na.rm = TRUE),
            Lower = quantile(Elasticity, 0.025, na.rm = TRUE),
            Upper = quantile(Elasticity, 0.975, na.rm = TRUE),
            .groups = "drop")


## Transient LTRE --------------------------------------------------------------

library(tidyverse)
library(data.table)
library(patchwork)
library(scales)

oneProp <- TRUE # whether age structure should be summed

LTREresults <- readRDS('results/LTREresults_random.rds')
plotFolder <- c("figures/resultsDave2Covs")
nAge <- 19

# extract relevant data
contData <- LTREresults$contData
nSamples <- length(LTREresults$contList$cont[[1]])

contData <- contData %>%
  mutate(draw = rep(1:(nrow(.) / length(unique(Variable))),
                    each = length(unique(Variable))))

if(oneProp){
  tmp <- contData %>% 
    filter(Variable %in% c("pYF", "pSA", "pAD_all")) %>% 
    group_by(draw) %>% 
    summarise(Contribution = sum(Contribution), .groups = "drop") %>% 
    mutate(Variable = "p_all", type = "Population structure")
  
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
  
  remove(tmp)
  
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

# pick colours
plotData$type <- factor(plotData$type, levels = types)
plot.colours <- paletteer::paletteer_c("grDevices::Temps", length(unique(plotData$type)))

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
        plot.margin = margin(1, 3, 1, 3)); c.sum

# ggsave("figures/resultsDave2Covs/LTREsum.jpeg", width = 20.0, height = 12.0, units = c("cm"), dpi = 600)

e.sum / c.sum

# ggsave("figures/resultsDave2Covs/elas&ltre.jpeg", width = 18.0, height = 18.0, units = c("cm"), dpi = 600)

# survival panel
S.colours <- c(plot.colours[3:4], rep(plot.colours[5], 8), rep(plot.colours[6], 10))
names(S.colours) <- c("sYF", "sSA", paste0("sAD_", 2:nAge))

p.S <- ggplot(subset(contData, Variable %in% c("sYF", "sSA", paste0("sAD_", 2:nAge)))) +
  geom_violin(aes(x = factor(Variable, levels = c("sYF", "sSA", paste0("sAD_", 2:nAge))),
                  y = Contribution, fill = Variable), alpha = 0.5, scale = "width", draw_quantiles = 0.5) +
  geom_hline(yintercept = 0, colour = "grey70", linetype = "dashed") +
  ylab("Contribution") +
  xlab("") +
  labs(title = "a) Survival") +
  scale_x_discrete(labels = expression(
    italic(S)[0], italic(S)[1],
    italic(S)[2], italic(S)[3], italic(S)[4], italic(S)[5], italic(S)[6],
    italic(S)[7], italic(S)[8], italic(S)[9], italic(S)[10], italic(S)[11],
    italic(S)[12], italic(S)[13], italic(S)[14], italic(S)[15], italic(S)[16],
    italic(S)[17], italic(S)[18], italic(S)[19]
  )) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  scale_fill_manual(values = S.colours) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.title = element_text(size = 10),
        plot.margin = margin(1, 3, 1, 3)); p.S

# birth rate panel
B.colours <- c(rep(plot.colours[1], 18))
names(B.colours) <- c(paste0("BR_", 2:nAge))

p.B <- ggplot(subset(contData, Variable %in% c(paste0("BR_", 2:nAge)))) +
  geom_violin(aes(x = factor(Variable, levels = c(paste0("BR_", 2:nAge))),
                  y = Contribution, fill = Variable), alpha = 0.5, scale = "width", draw_quantiles = 0.5) +
  geom_hline(yintercept = 0, colour = "grey70", linetype = "dashed") +
  ylab("Contribution") +
  xlab("") +
  labs(title = "b) Birth rate") +
  scale_x_discrete(labels = expression(
    italic(B)[2], italic(B)[3], italic(B)[4], italic(B)[5], italic(B)[6],
    italic(B)[7], italic(B)[8], italic(B)[9], italic(B)[10], italic(B)[11],
    italic(B)[12], italic(B)[13], italic(B)[14], italic(B)[15], italic(B)[16],
    italic(B)[17], italic(B)[18], italic(B)[19]
  )) +
  scale_y_continuous(labels = scales::label_number(),
                     expand = expansion(mult = c(0, 0.02))) +
  scale_fill_manual(values = B.colours) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.title = element_text(size = 10),
        plot.margin = margin(1, 3, 1, 3)); p.B

# reproductive success panel
R.colours <- c(rep(plot.colours[2], 18))
names(R.colours) <- c(paste0("sPY_", 2:nAge))

p.R <- ggplot(subset(contData, Variable %in% c(paste0("sPY_", 2:nAge)))) +
  geom_violin(aes(x = factor(Variable, levels = c(paste0("sPY_", 2:nAge))),
                  y = Contribution, fill = Variable), alpha = 0.5, scale = "width", draw_quantiles = 0.5) +
  geom_hline(yintercept = 0, colour = "grey70", linetype = "dashed") +
  ylab("Contribution") +
  xlab("") +
  labs(title = "c) Survival of pouch young") +
  scale_x_discrete(labels = expression(
    italic(SP)[2], italic(SP)[3], italic(SP)[4], italic(SP)[5], italic(SP)[6],
    italic(SP)[7], italic(SP)[8], italic(SP)[9], italic(SP)[10], italic(SP)[11],
    italic(SP)[12], italic(SP)[13], italic(SP)[14], italic(SP)[15], italic(SP)[16],
    italic(SP)[17], italic(SP)[18], italic(SP)[19]
  )) +
  scale_y_continuous(limits = c(0, 0.001), expand = expansion(mult = c(0, 0.02))) +
  scale_fill_manual(values = R.colours) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.title = element_text(size = 10),
        plot.margin = margin(1, 3, 1, 3)); p.R

# population structure panel
P.colours <- c(rep(plot.colours[7], 20))
names(P.colours) <- c("pYF", "pSA", paste0("pAD_", 2:nAge))

p.P <- ggplot(subset(contData, Variable %in% c("pYF", "pSA", paste0("pAD_", 2:nAge)))) +
  geom_violin(aes(x = factor(Variable, levels = c("pYF", "pSA", paste0("pAD_", 2:nAge))),
                  y = Contribution, fill = Variable), alpha = 0.5, scale = "width", draw_quantiles = 0.5) +
  geom_hline(yintercept = 0, colour = "grey70", linetype = "dashed") +
  ylab("Contribution") +
  xlab("") +
  labs(title = "d) Population proportions") +
  scale_x_discrete(labels = expression(
    italic(P)[0], italic(P)[1],
    italic(P)[2], italic(P)[3], italic(P)[4], italic(P)[5], italic(P)[6],
    italic(P)[7], italic(P)[8], italic(P)[9], italic(P)[10], italic(P)[11],
    italic(P)[12], italic(P)[13], italic(P)[14], italic(P)[15], italic(P)[16],
    italic(P)[17], italic(P)[18], italic(P)[19]
  )) +
  scale_y_continuous(limits = c(-0.01, 0.01), expand = expansion(mult = c(0, 0.02))) +
  scale_fill_manual(values = P.colours) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.title = element_text(size = 10),
        plot.margin = margin(1, 3, 1, 3)); p.P

# combine panels
# ((p.S + labs(tag = "a)")) / (p.B + labs(tag = "b)")) / (p.R + labs(tag = "c)")) / (p.P + labs(tag = "d)")))
p.S / p.B / p.R / p.P

# ggsave("figures/resultsDave2Covs/LTREage2.jpeg", width = 20.0, height = 24.0, units = c("cm"), dpi = 600)

# summaries to report
LTREsummary <- plotData %>%
  group_by(Variable) %>%
  summarise(Mean  = mean(Contribution, na.rm = TRUE),
            Lower = quantile(Contribution, 0.025, na.rm = TRUE),
            Upper = quantile(Contribution, 0.975, na.rm = TRUE),
            .groups = "drop")


## Density data ----------------------------------------------------------------

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

# ggsave("figures/densData_shrunkCIs.jpeg", width = 18.0, height = 10.0, units = c("cm"), dpi = 600)

tmp <- as.data.frame(cbind(dave$year, dave$dens, heloise$dens)) %>% 
  rename(year = 'V1', dave = 'V2', heloise = 'V3')

library(ggpubr)
ggscatter(tmp, x = "dave", y = "heloise",
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "pearson",
          alpha = 0.2, ggtheme = theme_bw())
