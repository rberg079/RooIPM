# 20 April 2026
# Plot results of IPM

## Set up ----------------------------------------------------------------------

library(tidyverse)
library(data.table)
library(patchwork)
library(scales)

nYear <- 17
nAge  <- 18

# load results
out.mcmc <- readRDS('results/IPM_CJSen_RSen_AB_DynDens_dCJS_12_noW_stochV.rds')
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
  geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "#C398B7", alpha = 0.4) +
  geom_line(color = "#673C5B", linewidth = 0.8) +
  scale_x_continuous(limits = c(2009, 2024),
                     breaks = c(2010, 2012, 2014, 2016, 2018, 2020, 2022, 2024)) +
  scale_y_continuous(breaks = pretty_breaks()) +
  labs(y = "Population size") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none"); pop

# ggsave("figures/results12ageCs/nTOT.jpeg", width = 18.0, height = 10.0, units = c("cm"), dpi = 600)


## Survival --------------------------------------------------------------------

# indices
S_idx  <- grep("^S\\[", colnames(out.mat))

# build summary dataframe
df <- expand.grid(Age = 1:13, Year = 1:16)
df$Mean <- apply(out.mat[, S_idx, drop = FALSE], 2, mean, na.rm = TRUE)
df$Lower <- apply(out.mat[, S_idx, drop = FALSE], 2, quantile, probs = 0.025, na.rm = TRUE)
df$Upper <- apply(out.mat[, S_idx, drop = FALSE], 2, quantile, probs = 0.975, na.rm = TRUE)

# select ages!
plotAges <- c(0, 1, 2, 4, 6, 8, 10, 12)

# pick colours
cols <- c(
  "0"  = "#F8756C",
  "1"  = "#CC9400",
  "2"  = "#7CAE00",
  "4"  = "#00BD65",
  "6"  = "#00BEC3",
  "8"  = "#00A8FF",
  "10" = "#C980FF",
  "12" = "#FF61CC"
)

# plot
surv <- df %>%
  mutate(Year = Year + 2007,
         Age  = factor(Age-1)) %>%
  filter(Age %in% plotAges) %>% 
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
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank()); surv

# ggsave("figures/results12ageCs/survival.jpeg", width = 18.0, height = 10.0, units = c("cm"), dpi = 600)


## Reproductive success --------------------------------------------------------

# indices
Bt_idx  <- grep("^Bt\\[", colnames(out.mat))
sPY_idx <- grep("^sPY\\[", colnames(out.mat))

# extract matrices
Bt  <- out.mat[, Bt_idx,  drop = FALSE]
Bt  <- Bt[, rep(1:ncol(Bt), each = 19)]
sPY <- out.mat[, sPY_idx, drop = FALSE]

# compute reproductive output
R <- 0.5 * Bt * sPY

# build summary dataframe
df <- expand.grid(Age = 1:19, Year = 1:16)
df$Mean <- apply(R, 2, mean, na.rm = TRUE)
df$Lower <- apply(R, 2, quantile, probs = 0.025, na.rm = TRUE)
df$Upper <- apply(R, 2, quantile, probs = 0.975, na.rm = TRUE)

# select ages!
plotAges <- c(2, 4, 6, 8, 10, 12)

# pick colours
cols <- c(
  "2"  = "#7CAE00",
  "4"  = "#00BD65",
  "6"  = "#00BEC3",
  "8"  = "#00A8FF",
  "10" = "#C980FF",
  "12" = "#FF61CC"
)

# plot
rout <- df %>%
  mutate(Year = Year + 2007,
         Age  = factor(Age)) %>%
  filter(Age %in% plotAges) %>% 
  ggplot(aes(x = Year, y = Mean, group = Age, colour = Age)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = Age), alpha = 0.2, colour = NA, show.legend = F) +
  geom_line(linewidth = 0.8, show.legend = F) +
  scale_colour_manual(values = cols) +
  scale_fill_manual(values = cols) +
  scale_x_continuous(limits = c(2008, 2024),
                     breaks = c(2008, 2010, 2012, 2014, 2016, 2018, 2020, 2022, 2024)) +
  scale_y_continuous(limits = c(0, 0.6), breaks = c(0.0, 0.2, 0.4, 0.6)) +
  labs(x = "Year", y = "Reproductive output", colour = "Age", fill = "Age") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none"); rout

# ggsave("figures/results12ageCs/routput.jpeg", width = 18.0, height = 10.0, units = c("cm"), dpi = 600)

# combine with survival plot
# surv / rout

# ggsave("figures/results12ageCs/surv&rout.jpeg", width = 18.0, height = 18.0, units = c("cm"), dpi = 600)

# ...& population size
(surv / rout / pop) +
  plot_layout(guides = "collect") +
  theme(legend.position = "right")

# ggsave("figures/results12ageCs/surv&rout&pop.jpeg", width = 18.0, height = 22.0, units = c("cm"), dpi = 600)


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
  make_df(bD, "Density"),
  make_df(bV, "Forage")
)

# plot
df %>%
  mutate(Age = factor(Age-1)) %>%
  ggplot(aes(x = x, y = Mean, colour = Age)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = Age), alpha = 0.2, colour = NA) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ covariate, scales = "free_x") +
  scale_y_continuous(limits = c(0.2, 1), breaks = c(0.2, 0.4, 0.6, 0.8, 1.0)) +
  labs(x = "Scaled covariate value", y = "Survival",
       colour = "Age", fill = "Age") +
  theme_bw()

# ggsave("figures/results12ageCs/coveffects.jpeg", width = 18.0, height = 10.0, units = c("cm"), dpi = 600)


## Covariate values ------------------------------------------------------------

# indices
D_idx  <- grep("^dens\\.true", colnames(out.mat))
V_idx  <- grep("^veg\\.true", colnames(out.mat))

# build summary dataframe
df <- expand.grid(Year = 1:17)

dens <- apply(out.mat[, D_idx, drop = FALSE], 2, mean, na.rm = TRUE)
veg  <- apply(out.mat[, V_idx, drop = FALSE], 2, mean, na.rm = TRUE)

dens <- scale(dens)[,1]

dens <- cbind(df, dens) %>% rename(value = dens) %>% mutate(covariate = "Density")
veg  <- cbind(df, veg) %>% rename(value = veg) %>% mutate(covariate = "Forage")

df <- rbind(dens, veg)

# plot
covs <- df %>%
  filter(Year > 1) %>% 
  mutate(Year = Year + 2007) %>%
  ggplot(aes(x = Year, y = value, colour = covariate)) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = 0, colour = "grey40") +
  scale_x_continuous(limits = c(2008, 2024),
                     breaks = c(2008, 2010, 2012, 2014, 2016, 2018, 2020, 2022, 2024)) +
  labs(y = "Scaled covariate value", colour = "Covariate") +
  theme_bw(); covs

# ggsave("figures/results12ageCs/covsVStime.jpeg", width = 18.0, height = 10.0, units = c("cm"), dpi = 600)

# combine with survival & population size
(surv / rout / covs) +
  plot_layout(guides = "collect") +
  theme(legend.position = "right")

# ggsave("figures/results12ageCs/surv&rout&covs.jpeg", width = 18.0, height = 22.0, units = c("cm"), dpi = 600)


## Lambda ----------------------------------------------------------------------

# load results
paramSamples <- readRDS('results/paramSamples.rds')

# extract lambda
lambda_mat <- paramSamples$t$lambda
lambda_vec <- as.vector(lambda_mat)

# summary stats
mean(lambda_vec)
quantile(lambda_vec, c(0.025, 0.5, 0.975))

# build summary dataframe
df <- data.frame(Year  = 1:16,
                 Mean  = apply(lambda_mat, 2, mean),
                 Lower = apply(lambda_mat, 2, quantile, 0.025),
                 Upper = apply(lambda_mat, 2, quantile, 0.975))

lambda <- df %>%
  mutate(Year = Year + 2007) %>%
  ggplot(aes(x = Year, y = Mean)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.4, fill = "#B6967C") +
  geom_line(colour = "#684F3B", linewidth = 0.8) +
  geom_hline(yintercept = 1, colour = "grey40") +
  scale_x_continuous(limits = c(2009, 2024),
                     breaks = c(2010, 2012, 2014, 2016, 2018, 2020, 2022, 2024)) +
  labs(x = "Year", y = expression("Population growth" ~ (lambda))) +
  theme_bw(); lambda

# ggsave("figures/results12ageCs/lambda.jpeg", width = 18.0, height = 10.0, units = c("cm"), dpi = 600)

# combine with pop size
pop / lambda
pAGE / lambda

# ggsave("figures/results12ageCs/pop&lambda.jpeg", width = 18.0, height = 18.0, units = c("cm"), dpi = 600)


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
cols <- c("0"   = "#E8D6CB",
          "1"   = "#D0ADA7",
          setNames(colorRampPalette(c("#AD6A6C", "#5D2E46"))(10), as.character(2:11)),
          "12+" = "#371B29")

# cols <- c("0"   = "#C9CBA3",
#           "1"   = "#FFE1A8",
#           setNames(colorRampPalette(c("#F1A782", "#8E494C"))(10), as.character(2:11)),
#           "12+" = "#723D46")

# bar plot
pAGE <- df %>%
  ggplot(aes(x = Year, y = N, fill = AgeGroup)) +
  # ggplot(aes(x = Year, y = Prop, fill = AgeGroup)) +
  geom_col(width = 0.8) +
  scale_fill_manual(values = cols) +
  scale_x_continuous(limits = c(2009, 2024),
                     breaks = c(2010, 2012, 2014, 2016, 2018, 2020, 2022, 2024)) +
  labs(x = "Year", y = "Proportion of the population", fill = "Age") +
  theme_bw(); pAGE

# ggsave("figures/results12ageCs/NsBars.jpeg", width = 18.0, height = 10.0, units = c("cm"), dpi = 600)

# ribbon plot
pAGE <- df %>% 
  arrange(Year, AgeGroup) %>%
  group_by(Year) %>%
  # mutate(ymin = cumsum(lag(N, default = 0)),
  #        ymax = cumsum(N)) %>%
  mutate(ymin = cumsum(lag(Prop, default = 0)),
         ymax = cumsum(Prop)) %>%
  ungroup() %>%
  ggplot(aes(x = Year, fill = AgeGroup)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), colour = NA) +
  scale_fill_manual(values = cols) +
  scale_x_continuous(limits = c(2009, 2024),
                     breaks = c(2010, 2012, 2014, 2016, 2018, 2020, 2022, 2024)) +
  # scale_y_continuous(limits = c(0, 750)) +
  labs(x = "Year", y = "Population size", fill = "Age") +
  theme_bw(); pAGE

# ggsave("figures/results12ageCs/NsRibbons.jpeg", width = 18.0, height = 10.0, units = c("cm"), dpi = 600)

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

