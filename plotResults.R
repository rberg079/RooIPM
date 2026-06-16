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
out.mcmc <- readRDS('results/IPM_CJSen_RSen_AB_DynDens_dCJS_12_noW_stochV_long_BR.rds')
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
  filter(Year > 1) %>% 
  mutate(Year = Year + 2007) %>% 
  ggplot(aes(x = Year, y = Mean)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "#62414B", alpha = 0.2) + #C398B7
  geom_line(color = "#62414B", linewidth = 0.8) + #673C5B
  scale_x_continuous(limits = c(2009, 2025),
                     breaks = c(2009, 2011, 2013, 2015, 2017, 2019, 2021, 2023, 2025)) +
  # scale_x_continuous(limits = c(2008, 2025),
  #                    breaks = c(2008, 2010, 2012, 2014, 2016, 2018, 2020, 2022, 2024)) +
  scale_y_continuous(limits = c(200, 800),
                     breaks = c(200, 400, 600, 800)) +
  labs(y = "Population size") +
  theme_bw() +
  theme(
        axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none"
        ); pop

# ggsave("figures/results25&BR/nTOT.jpeg", width = 18.0, height = 10.0, units = c("cm"), dpi = 600)


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
pop <- df %>% 
  filter(Year > 1) %>% 
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
  scale_x_continuous(limits = c(2009, 2025),
                     breaks = c(2009, 2011, 2013, 2015, 2017, 2019, 2021, 2023, 2025)) +
  scale_y_continuous(limits = c(0, 800), breaks = scales::pretty_breaks()) +
  labs(y = "Population size", colour = "Age group", fill = "Age group") +
  theme_bw(); pop

# ggsave("figures/results25&BR/allNs.jpeg", width = 18.0, height = 10.0, units = c("cm"), dpi = 600)


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
#   "0"  = "#F8756C",
#   "1"  = "#CC9400",
#   "2"  = "#7CAE00",
#   "4"  = "#00BD65",
#   "6"  = "#00BEC3",
#   "8"  = "#00A8FF",
#   "10" = "#C980FF",
#   "12+" = "#FF61CC"
# )
# 
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

cols <- c(
  "0"  = "#8E5EA2",
  "1"  = "#277DA1",
  "2"  = "#00A896",
  "4"  = "#F9C74F",
  "6"  = "#F8A23A",
  "8"  = "#F3722C",
  "10" = "#D62828",
  "12+" = "#D96C9D"
)

# cols <- c(
#   "0"  = "#283D3B",
#   "1"  = "#197278",
#   "2"  = "#83A8A6",
#   "4"  = "#EDDDD4",
#   "6"  = "#AE9D96",
#   "8"  = "#D99185",
#   "10" = "#C44536",
#   "12+" = "#772E25"
# )

# plot
surv <- df %>%
  mutate(Year = Year + 2007,
         Age  = factor(Age-1,
                       levels = c(0, 1, 2, 4, 6, 8, 10, 12),
                       labels = c("0", "1", "2", "4", "6", "8", "10", "12+"))) %>%
  filter(!is.na(Age)) %>% 
  ggplot(aes(x = Year, y = Mean, group = Age, colour = Age)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = Age), alpha = 0.2, colour = NA) +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = cols) +
  scale_fill_manual(values = cols) +
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

# ggsave("figures/results25&BR/survival.jpeg", width = 18.0, height = 10.0, units = c("cm"), dpi = 600)


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

# pick colours
# cols <- c(
#   "2"  = "#7CAE00",
#   "4"  = "#00BD65",
#   "6"  = "#00BEC3",
#   "8"  = "#00A8FF",
#   "10" = "#C980FF",
#   "12+" = "#FF61CC"
# )
# 
# cols <- c(
#   "2"  = "#DE612B",
#   "4"  = "#F4A261",
#   "6"  = "#E9C46A",
#   "8"  = "#2A9D8F",
#   "10" = "#264653",
#   "12+" = "#936271"
# )

cols <- c(
  "2"  = "#00A896",
  "4"  = "#F9C74F",
  "6"  = "#F8A23A",
  "8"  = "#F3722C",
  "10" = "#D62828",
  "12+" = "#D96C9D"
)

# plot
rs <- df %>%
  mutate(Year = Year + 2007,
         Age  = factor(Age,
                       levels = c(2, 4, 6, 8, 10, 12),
                       labels = c("2", "4", "6", "8", "10", "12+"))) %>%
  filter(!is.na(Age)) %>% 
  ggplot(aes(x = Year, y = Mean, group = Age, colour = Age)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = Age), alpha = 0.2, colour = NA, show.legend = F) +
  geom_line(linewidth = 0.8, show.legend = F) +
  scale_colour_manual(values = cols) +
  scale_fill_manual(values = cols) +
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
    # axis.ticks.x = element_blank()
  ); rs

# ggsave("figures/results25&BR/PYsurv.jpeg", width = 18.0, height = 10.0, units = c("cm"), dpi = 600)

# combine with survival plot
# surv / rs

# combine birth rate & PY survival
(br / rs) +
  plot_layout(guides = "collect") +
  theme(legend.position = "right")

# ggsave("figures/results25&BR/birth&PYsurv.jpeg", width = 18.0, height = 18.0, units = c("cm"), dpi = 600)

# ...& population size
(surv / rs / pop) +
  plot_layout(guides = "collect") +
  theme(legend.position = "right")

# ggsave("figures/results25&BR/surv&rs&pop.jpeg", width = 18.0, height = 22.0, units = c("cm"), dpi = 600)


## Covariate effects -----------------------------------------------------------

inv_logit <- function(x) 1 / (1 + exp(-x))

keep_ages <- c(1, 11, 12, 13)

# extract betas
bD <- out.mat[, grep("^BetaD\\.S", colnames(out.mat))]
bV <- out.mat[, grep("^BetaV\\.S", colnames(out.mat))]

# extract baseline intercepts
mus  <- paste0("Mu.S[", keep_ages, "]")
base <- out.mat[, mus, drop = FALSE]

# covariate sequence
x <- seq(-2, 2, length.out = 50)

# build summary dfs
make_df <- function(betas, covariate) {
  results <- list()
  
  for (i in seq_along(keep_ages)) {
    age <- keep_ages[i]
    intercepts <- base[, i] 
    
    for (xx in x) {
      preds <- inv_logit(intercepts + betas * xx)
      results[[length(results) + 1]] <- data.frame(
        Age       = age,
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

df <- rbind(
  make_df(bD, "Population density"),
  make_df(bV, "Forage availability")
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

cols <- c(
  "0"   = "#7C528E",
  "10"  = "#BE2323",
  "11"  = "#696969",
  "12+" = "#D96C9D"
)

# plot
df %>%
  mutate(Age = factor(Age - 1,
                      levels = c(0, 10, 11, 12),
                      labels = c("0", "10", "11", "12+"))) %>%
  ggplot(aes(x = x, y = Mean, colour = Age)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = Age), alpha = 0.2, colour = NA) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ covariate, scales = "free_x") +
  scale_colour_manual(values = cols) +
  scale_fill_manual(values = cols) +
  scale_y_continuous(limits = c(0.2, 1), breaks = c(0.2, 0.4, 0.6, 0.8, 1.0)) +
  labs(x = "Scaled covariate value", y = "Survival",
       colour = "Age", fill = "Age") +
  theme_bw()

# ggsave("figures/results25&BR/coveffects.jpeg", width = 18.0, height = 10.0, units = c("cm"), dpi = 600)


## Covariate values ------------------------------------------------------------

# indices
D_idx  <- grep("^dens\\.true", colnames(out.mat))[1:17]
V_idx  <- grep("^veg\\.true", colnames(out.mat))

# build summary dataframe
df <- expand.grid(Year = 1:17)

dens <- apply(out.mat[, D_idx, drop = FALSE], 2, mean, na.rm = TRUE)
veg  <- apply(out.mat[, V_idx, drop = FALSE], 2, mean, na.rm = TRUE)

dens <- scale(dens)[,1]

dens <- cbind(df, dens) %>% rename(value = dens) %>% mutate(covariate = "Population density")
veg  <- cbind(df, veg) %>% rename(value = veg) %>% mutate(covariate = "Forage availability")

df <- rbind(dens, veg)

# pick colours
cols <- c(
  "Population density" = "#335B5B",
  "Forage availability"  = "#75A366"
)

# plot
covs <- df %>%
  filter(Year > 1) %>% 
  mutate(Year = Year + 2007) %>%
  ggplot(aes(x = Year, y = value, colour = covariate)) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = 0, colour = "grey40") +
  scale_colour_manual(values = cols) +
  scale_x_continuous(limits = c(2008, 2024),
                     breaks = c(2008, 2010, 2012, 2014, 2016, 2018, 2020, 2022, 2024)) +
  labs(y = "Scaled covariate value", colour = "Covariate") +
  theme_bw(); covs

# ggsave("figures/results25&BR/covsVStime.jpeg", width = 18.0, height = 10.0, units = c("cm"), dpi = 600)

# combine with survival & reproductive output
(surv / rs / covs) +
  plot_layout(guides = "collect") +
  theme(legend.position = "right")

# ggsave("figures/results25&BR/surv&rs&covs.jpeg", width = 20.0, height = 22.0, units = c("cm"), dpi = 600)


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
  # scale_x_continuous(limits = c(2008, 2025),
  #                    breaks = c(2008, 2010, 2012, 2014, 2016, 2018, 2020, 2022, 2024)) +
  scale_x_continuous(limits = c(2009, 2025),
                     breaks = c(2009, 2011, 2013, 2015, 2017, 2019, 2021, 2023, 2025)) +
  scale_y_continuous(limits = c(NA, 1.4),
                     breaks = c(0.6, 0.8, 1.0, 1.2, 1.4)) +
  labs(x = "Year", y = expression("Population growth rate" ~ (lambda))) +
  theme_bw(); lambda

# ggsave("figures/results25&BR/lambda.jpeg", width = 18.0, height = 10.0, units = c("cm"), dpi = 600)

# combine with pop size
pop / lambda
pAGE / lambda

# ggsave("figures/results25&BR/pop&lambda.jpeg", width = 18.0, height = 18.0, units = c("cm"), dpi = 600)


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

cols <- c("0"   = "#C9CBA3",
          "1"   = "#FFE1A8",
          setNames(colorRampPalette(c("#EDA297", "#DD5540"))(8), as.character(2:9)),
          setNames(colorRampPalette(c("#854752", "#5D3239"))(3), as.character(c("10", "11", "12+"))))

cols <- c("0"   = "#E3AD78",
          "1"   = "#EFD6AC",
          setNames(colorRampPalette(c("#A0A290", "#4D4F40"))(8), as.character(2:9)),
          setNames(colorRampPalette(c("#695C70", "#4A404F"))(3), as.character(c("10", "11", "12+"))))

cols <- c("0"   = "#CA562C",
          "1"   = "#DE8A5A",
          setNames(colorRampPalette(c("#FAF5DB", "#F2D4A4"))(3), as.character(2:4)),
          setNames(colorRampPalette(c("#B5B991", "#657359"))(5), as.character(5:9)),
          setNames(colorRampPalette(c("#324935", "#19241A"))(3), as.character(c(10:11, "12+"))))

cols <- c("0"   = "#671313",
          "1"   = "#752915",
          setNames(colorRampPalette(c("#9E691A", "#E7D18B"))(3), as.character(2:4)),
          setNames(colorRampPalette(c("#FFF3B0", "#798C69"))(5), as.character(5:9)),
          setNames(colorRampPalette(c("#426061", "#1C3236"))(3), as.character(c(10:11, "12+"))))

# cols <- paletteer::paletteer_c("grDevices::Temps", 7)

# ribbon plot
pAGE <- df %>% 
  arrange(Year, desc(AgeGroup)) %>% # desc() to flip order
  group_by(Year) %>%
  # mutate(ymin = cumsum(lag(Prop, default = 0)),
  #        ymax = cumsum(Prop)) %>%
  mutate(ymin = cumsum(lag(N, default = 0)),
  ymax = cumsum(N)) %>%
  ungroup() %>%
  ggplot(aes(x = Year, fill = AgeGroup)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), colour = NA) +
  scale_fill_manual(values = cols) +
  # guides(fill = guide_legend(reverse = TRUE)) +
  scale_x_continuous(limits = c(2009, 2025),
                     breaks = c(seq(2009, 2025, by = 2))) +
  scale_y_continuous(limits = c(0, 750)) +
  labs(x = "Year", y = "Population size", fill = "Age") +
  theme_bw(); pAGE

# ggsave("figures/results25&BR/NsRibbons2.jpeg", width = 18.0, height = 10.0, units = c("cm"), dpi = 600)

# bar plot
pAGE <- df %>%
  ggplot(aes(x = Year, y = N, fill = AgeGroup)) +
  geom_col(width = 0.8) +
  scale_fill_manual(values = cols) +
  scale_x_continuous(limits = c(2009, 2024),
                     breaks = c(2010, 2012, 2014, 2016, 2018, 2020, 2022, 2024)) +
  labs(x = "Year", y = "Proportion of the population", fill = "Age") +
  theme_bw(); pAGE

# ggsave("figures/results25&BR/NsBars.jpeg", width = 18.0, height = 10.0, units = c("cm"), dpi = 600)

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

# ggsave("figures/results25&BR/SENSsum.jpeg", width = 20.0, height = 12.0, units = c("cm"), dpi = 600)

# birth rate panel
B.colours <- c(rep(plot.colours[1], 18))
names(B.colours) <- c(paste0("BR_", 2:nAge))

p.B <- ggplot(subset(sensData, Variable %in% c(paste0("BR_", 2:nAge)))) +
  geom_violin(aes(x = factor(Variable, levels = c(paste0("BR_", 2:nAge))),
                  y = Sensitivity, fill = Variable), alpha = 0.5, scale = "width", draw_quantiles = 0.5) +
  geom_hline(yintercept = 0, colour = "grey70", linetype = "dashed") +
  ylab("Sensitivity") +
  xlab("") +
  scale_x_discrete(labels = expression(B[2], B[3], B[4], B[5], B[6], B[7], B[8], B[9], B[10],
                                       B[11], B[12], B[13], B[14], B[15], B[16], B[17], B[18])) +
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
  scale_x_discrete(labels = expression(SP[2], SP[3], SP[4], SP[5], SP[6], SP[7], SP[8], SP[9], SP[10],
                                       SP[11], SP[12], SP[13], SP[14], SP[15], SP[16], SP[17], SP[18])) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  scale_fill_manual(values = R.colours) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.title = element_text(size = 10),
        plot.margin = margin(1, 3, 1, 3)); p.R

# survival panel
S.colours <- c(plot.colours[3:4], rep(plot.colours[5], 8), rep(plot.colours[6], 10))
names(S.colours) <- c("sYF", "sSA", paste0("sAD_", 2:nAge))

p.S <- ggplot(subset(sensData, Variable %in% c("sYF", "sSA", paste0("sAD_", 2:nAge)))) +
  geom_violin(aes(x = factor(Variable, levels = c("sYF", "sSA", paste0("sAD_", 2:nAge))),
                  y = Sensitivity, fill = Variable), alpha = 0.5, scale = "width", draw_quantiles = 0.5) +
  geom_hline(yintercept = 0, colour = "grey70", linetype = "dashed") +
  ylab("Sensitivity") +
  xlab("") +
  scale_x_discrete(labels = expression(S[0], S[1],
                                       S[2], S[3], S[4], S[5], S[6], 
                                       S[7], S[8], S[9], S[10], S[11], S[12],
                                       S[13], S[14], S[15], S[16], S[17], S[18], S[19])) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  scale_fill_manual(values = S.colours) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.title = element_text(size = 10),
        plot.margin = margin(1, 3, 1, 3)); p.S

# population structure panel
P.colours <- c(rep(plot.colours[7], 20))
names(P.colours) <- c("pYF", "pSA", paste0("pAD_", 2:nAge))

p.P <- ggplot(subset(sensData, Variable %in% c("pYF", "pSA", paste0("pAD_", 2:nAge)))) +
  geom_violin(aes(x = factor(Variable, levels = c("pYF", "pSA", paste0("pAD_", 2:nAge))),
                  y = Sensitivity, fill = Variable), alpha = 0.5, scale = "width", draw_quantiles = 0.5) +
  geom_hline(yintercept = 0, colour = "grey70", linetype = "dashed") +
  ylab("Sensitivity") +
  xlab("") +
  scale_x_discrete(labels = expression(P[0], P[1],
                                       P[2], P[3], P[4], P[5], P[6], 
                                       P[7], P[8], P[9], P[10], P[11], P[12],
                                       P[13], P[14], P[15], P[16], P[17], P[18], P[19])) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  scale_fill_manual(values = P.colours) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.title = element_text(size = 10),
        plot.margin = margin(1, 3, 1, 3)); p.P

# combine panels
(p.sum + labs(tag = "a)")) / ((p.B + labs(tag = "b)")) / (p.R + labs(tag = "c)")) / (p.S + labs(tag = "d)")) / (p.P + labs(tag = "e)"))) +
  plot_layout(heights = c(0.4, 0.6))

((p.S + labs(tag = "a)")) / (p.B + labs(tag = "b)")) / (p.R + labs(tag = "c)")) / (p.P + labs(tag = "d)")))

# ggsave("figures/results25&BR/SENSage.jpeg", width = 20.0, height = 24.0, units = c("cm"), dpi = 600)

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

# ggsave("figures/results25&BR/ELASsum.jpeg", width = 20.0, height = 12.0, units = c("cm"), dpi = 600)

# birth rate panel
B.colours <- c(rep(plot.colours[1], 18))
names(B.colours) <- c(paste0("BR_", 2:nAge))

p.B <- ggplot(subset(elasData, Variable %in% c(paste0("BR_", 2:nAge)))) +
  geom_violin(aes(x = factor(Variable, levels = c(paste0("BR_", 2:nAge))),
                  y = Elasticity, fill = Variable), alpha = 0.5, scale = "width", draw_quantiles = 0.5) +
  geom_hline(yintercept = 0, colour = "grey70", linetype = "dashed") +
  ylab("Elasticity") +
  xlab("") +
  scale_x_discrete(labels = expression(B[2], B[3], B[4], B[5], B[6], B[7], B[8], B[9], B[10],
                                       B[11], B[12], B[13], B[14], B[15], B[16], B[17], B[18])) +
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

p.R <- ggplot(subset(elasData, Variable %in% c(paste0("sPY_", 2:nAge)))) +
  geom_violin(aes(x = factor(Variable, levels = c(paste0("sPY_", 2:nAge))),
                  y = Elasticity, fill = Variable), alpha = 0.5, scale = "width", draw_quantiles = 0.5) +
  geom_hline(yintercept = 0, colour = "grey70", linetype = "dashed") +
  ylab("Elasticity") +
  xlab("") +
  scale_x_discrete(labels = expression(SP[2], SP[3], SP[4], SP[5], SP[6], SP[7], SP[8], SP[9], SP[10],
                                       SP[11], SP[12], SP[13], SP[14], SP[15], SP[16], SP[17], SP[18])) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  scale_fill_manual(values = R.colours) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.title = element_text(size = 10),
        plot.margin = margin(1, 3, 1, 3)); p.R

# survival panel
S.colours <- c(plot.colours[3:4], rep(plot.colours[5], 8), rep(plot.colours[6], 10))
names(S.colours) <- c("sYF", "sSA", paste0("sAD_", 2:nAge))

p.S <- ggplot(subset(elasData, Variable %in% c("sYF", "sSA", paste0("sAD_", 2:nAge)))) +
  geom_violin(aes(x = factor(Variable, levels = c("sYF", "sSA", paste0("sAD_", 2:nAge))),
                  y = Elasticity, fill = Variable), alpha = 0.5, scale = "width", draw_quantiles = 0.5) +
  geom_hline(yintercept = 0, colour = "grey70", linetype = "dashed") +
  ylab("Elasticity") +
  xlab("") +
  scale_x_discrete(labels = expression(S[0], S[1],
                                       S[2], S[3], S[4], S[5], S[6], 
                                       S[7], S[8], S[9], S[10], S[11], S[12],
                                       S[13], S[14], S[15], S[16], S[17], S[18], S[19])) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  scale_fill_manual(values = S.colours) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.title = element_text(size = 10),
        plot.margin = margin(1, 3, 1, 3)); p.S

# population structure panel
P.colours <- c(rep(plot.colours[7], 20))
names(P.colours) <- c("pYF", "pSA", paste0("pAD_", 2:nAge))

p.P <- ggplot(subset(elasData, Variable %in% c("pYF", "pSA", paste0("pAD_", 2:nAge)))) +
  geom_violin(aes(x = factor(Variable, levels = c("pYF", "pSA", paste0("pAD_", 2:nAge))),
                  y = Elasticity, fill = Variable), alpha = 0.5, scale = "width", draw_quantiles = 0.5) +
  geom_hline(yintercept = 0, colour = "grey70", linetype = "dashed") +
  ylab("Elasticity") +
  xlab("") +
  scale_x_discrete(labels = expression(P[0], P[1],
                                       P[2], P[3], P[4], P[5], P[6], 
                                       P[7], P[8], P[9], P[10], P[11], P[12],
                                       P[13], P[14], P[15], P[16], P[17], P[18], P[19])) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  scale_fill_manual(values = P.colours) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.title = element_text(size = 10),
        plot.margin = margin(1, 3, 1, 3)); p.P

# combine panels
(p.sum + labs(tag = "a)")) / ((p.B + labs(tag = "b)")) / (p.R + labs(tag = "c)")) / (p.S + labs(tag = "d)")) / (p.P + labs(tag = "e)"))) +
  plot_layout(heights = c(0.4, 0.6))

((p.S + labs(tag = "a)")) / (p.B + labs(tag = "b)")) / (p.R + labs(tag = "c)")) / (p.P + labs(tag = "d)")))

# ggsave("figures/results25&BR/ELASage.jpeg", width = 20.0, height = 24.0, units = c("cm"), dpi = 600)

# summaries to report
elasSummary <- elasData %>%
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
plotFolder <- c("figures/results25&BR")
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
  scale_y_continuous(limits = c(-0.005, 0.008),
                     expand = expansion(mult = c(0, 0.02)),
                     breaks = c(-0.006, -0.004, -0.002, 0.00, 0.002, 0.004, 0.006, 0.008)) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.title = element_text(size = 10),
        plot.margin = margin(1, 3, 1, 3)); c.sum

# ggsave("figures/results25&BR/LTREsum.jpeg", width = 20.0, height = 12.0, units = c("cm"), dpi = 600)

e.sum / c.sum

# ggsave("figures/results25&BR/elas&ltre.jpeg", width = 18.0, height = 18.0, units = c("cm"), dpi = 600)

# birth rate panel
B.colours <- c(rep(plot.colours[1], 18))
names(B.colours) <- c(paste0("BR_", 2:nAge))

p.B <- ggplot(subset(contData, Variable %in% c(paste0("BR_", 2:nAge)))) +
  geom_violin(aes(x = factor(Variable, levels = c(paste0("BR_", 2:nAge))),
                  y = Contribution, fill = Variable), alpha = 0.5, scale = "width", draw_quantiles = 0.5) +
  geom_hline(yintercept = 0, colour = "grey70", linetype = "dashed") +
  ylab("Contribution") +
  xlab("") +
  scale_x_discrete(labels = expression(B[2], B[3], B[4], B[5], B[6], B[7], B[8], B[9], B[10],
                                       B[11], B[12], B[13], B[14], B[15], B[16], B[17], B[18])) +
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

p.R <- ggplot(subset(contData, Variable %in% c(paste0("sPY_", 2:nAge)))) +
  geom_violin(aes(x = factor(Variable, levels = c(paste0("sPY_", 2:nAge))),
                  y = Contribution, fill = Variable), alpha = 0.5, scale = "width", draw_quantiles = 0.5) +
  geom_hline(yintercept = 0, colour = "grey70", linetype = "dashed") +
  ylab("Contribution") +
  xlab("") +
  scale_x_discrete(labels = expression(SP[2], SP[3], SP[4], SP[5], SP[6], SP[7], SP[8], SP[9], SP[10],
                                       SP[11], SP[12], SP[13], SP[14], SP[15], SP[16], SP[17], SP[18])) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  scale_fill_manual(values = R.colours) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.title = element_text(size = 10),
        plot.margin = margin(1, 3, 1, 3)); p.R

# survival panel
S.colours <- c(plot.colours[3:4], rep(plot.colours[5], 8), rep(plot.colours[6], 10))
names(S.colours) <- c("sYF", "sSA", paste0("sAD_", 2:nAge))

p.S <- ggplot(subset(contData, Variable %in% c("sYF", "sSA", paste0("sAD_", 2:nAge)))) +
  geom_violin(aes(x = factor(Variable, levels = c("sYF", "sSA", paste0("sAD_", 2:nAge))),
                  y = Contribution, fill = Variable), alpha = 0.5, scale = "width", draw_quantiles = 0.5) +
  geom_hline(yintercept = 0, colour = "grey70", linetype = "dashed") +
  ylab("Contribution") +
  xlab("") +
  scale_x_discrete(labels = expression(S[0], S[1],
                                       S[2], S[3], S[4], S[5], S[6], 
                                       S[7], S[8], S[9], S[10], S[11], S[12],
                                       S[13], S[14], S[15], S[16], S[17], S[18], S[19])) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  scale_fill_manual(values = S.colours) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.title = element_text(size = 10),
        plot.margin = margin(1, 3, 1, 3)); p.S

# population structure panel
P.colours <- c(rep(plot.colours[7], 20))
names(P.colours) <- c("pYF", "pSA", paste0("pAD_", 2:nAge))

p.P <- ggplot(subset(contData, Variable %in% c("pYF", "pSA", paste0("pAD_", 2:nAge)))) +
  geom_violin(aes(x = factor(Variable, levels = c("pYF", "pSA", paste0("pAD_", 2:nAge))),
                  y = Contribution, fill = Variable), alpha = 0.5, scale = "width", draw_quantiles = 0.5) +
  geom_hline(yintercept = 0, colour = "grey70", linetype = "dashed") +
  ylab("Contribution") +
  xlab("") +
  scale_x_discrete(labels = expression(P[0], P[1],
                                       P[2], P[3], P[4], P[5], P[6], 
                                       P[7], P[8], P[9], P[10], P[11], P[12],
                                       P[13], P[14], P[15], P[16], P[17], P[18], P[19])) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  scale_fill_manual(values = P.colours) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.title = element_text(size = 10),
        plot.margin = margin(1, 3, 1, 3)); p.P

# combine panels
(p.sum + labs(tag = "a)")) / ((p.B + labs(tag = "b)")) / (p.R + labs(tag = "c)")) / (p.S + labs(tag = "d)")) / (p.P + labs(tag = "e)"))) +
  plot_layout(heights = c(0.4, 0.6))

((p.S + labs(tag = "a)")) / (p.B + labs(tag = "b)")) / (p.R + labs(tag = "c)")) / (p.P + labs(tag = "d)")))

# ggsave("figures/results25&BR/LTREage.jpeg", width = 20.0, height = 24.0, units = c("cm"), dpi = 600)

# summaries to report
LTREsummary <- contData %>%
  group_by(Variable) %>%
  summarise(Mean  = mean(Contribution, na.rm = TRUE),
            Lower = quantile(Contribution, 0.025, na.rm = TRUE),
            Upper = quantile(Contribution, 0.975, na.rm = TRUE),
            .groups = "drop")


## Density data ----------------------------------------------------------------

# wrangle Dave's data
source('R/wrangleData_en.R')
enData <- wrangleData_en(dens.data = "data/WPNP_Methods_Results_January2026.xlsx",
                         veg.data  = "data/biomass data April 2009 - July 2025_updated Feb2026.xlsx",
                         wea.data  = "data/Prom_Weather_2008-2023_updated Jan2026 RB.xlsx",
                         wind.data = "data/POWER_Point_Daily_20080101_20260331_10M.csv",
                         obs.data  = "data/PromObs_2008-2024.xlsx",
                         list      = "data/PromlistAllNov25.xlsx",
                         Dave      = TRUE)

year <- seq(2008, 2025, 1)
dave <- as.data.frame(cbind(year, enData$dens, enData$densE)) %>% 
  rename(dens = 'V2', densE = 'V3') %>% 
  mutate(data = "Dave") %>% 
  filter(year > 2008)

# wrangle Heloise's data
enData <- wrangleData_en(dens.data = "data/abundanceData_Proteus.csv",
                         veg.data  = "data/biomass data April 2009 - July 2025_updated Feb2026.xlsx",
                         wea.data  = "data/Prom_Weather_2008-2023_updated Jan2026 RB.xlsx",
                         wind.data = "data/POWER_Point_Daily_20080101_20260331_10M.csv",
                         obs.data  = "data/PromObs_2008-2024.xlsx",
                         list      = "data/PromlistAllNov25.xlsx",
                         Dave      = FALSE)

heloise <- as.data.frame(cbind(year, enData$dens, enData$densE)) %>% 
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

# ggsave("figures/densData.jpeg", width = 18.0, height = 10.0, units = c("cm"), dpi = 600)

