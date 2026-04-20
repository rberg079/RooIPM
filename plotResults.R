# 20 April 2026
# Plot results of IPM


## Set up ----------------------------------------------------------------------

library(tidyverse)
library(data.table)
library(scales)

nYear <- 17
nAge  <- 18

# load results
out.mcmc <- readRDS('results/IPM_CJSen_RSen_AB_DynDens_dCJS_12.rds')
out.mat <- do.call(rbind, lapply(out.mcmc, as.matrix))


## Plot population model -------------------------------------------------------

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
df %>% 
  filter(Year > 1) %>% 
  mutate(Year = Year + 2007) %>% 
  ggplot(aes(x = Year, y = Mean)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "#C398B7", alpha = 0.4) +
  geom_line(color = "#673C5B", linewidth = 0.8) +
  scale_x_continuous(breaks = c(2008, 2010, 2012, 2014, 2016, 2018, 2020, 2022, 2024)) +
  scale_y_continuous(breaks = pretty_breaks()) +
  labs(y = "Female population size") +
  theme_bw()

# ggsave("figures/results12ageCs/nTOT.jpeg", width = 18.0, height = 10.0, units = c("cm"), dpi = 600)


## Plot survival model ---------------------------------------------------------

# indices
S_idx  <- grep("^S\\[", colnames(out.mat))

# build summary dataframe
df <- expand.grid(Age = 1:13, Year = 1:16)
df$Mean <- apply(out.mat[, S_idx, drop = FALSE], 2, mean, na.rm = TRUE)
df$Lower <- apply(out.mat[, S_idx, drop = FALSE], 2, quantile, probs = 0.025, na.rm = TRUE)
df$Upper <- apply(out.mat[, S_idx, drop = FALSE], 2, quantile, probs = 0.975, na.rm = TRUE)

# select ages!
plotAges <- c(0, 1, 2, 4, 6, 8, 10, 12)

# plot
df %>%
  mutate(Year = Year + 2007,
         Age  = factor(Age-1)) %>%
  filter(Age %in% plotAges) %>% 
  ggplot(aes(x = Year, y = Mean, group = Age, colour = Age)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = Age), alpha = 0.2, colour = NA) +
  geom_line(linewidth = 0.8) +
  scale_x_continuous(breaks = c(2008, 2010, 2012, 2014, 2016, 2018, 2020, 2022)) +
  scale_y_continuous(limits = c(0,1), breaks = pretty_breaks()) +
  labs(x = "Year", y = "Survival", colour = "Age", fill = "Age") +
  theme_bw()

# ggsave("figures/results12ageCs/survival.jpeg", width = 18.0, height = 10.0, units = c("cm"), dpi = 600)


## Plots reproductive output ---------------------------------------------------

# indices
Bt_idx  <- grep("^Bt\\[", colnames(out.mat))
sPY_idx <- grep("^sPY\\[", colnames(out.mat))

# extract matrices
Bt  <- out.mat[, Bt_idx,  drop = FALSE]
Bt  <- Bt[, rep(1:ncol(Bt), each = 18)]
sPY <- out.mat[, sPY_idx, drop = FALSE]

# compute reproductive output
R <- 0.5 * Bt * sPY

# build summary dataframe
df <- expand.grid(Age = 1:18, Year = 1:16)
df$Mean <- apply(RO, 2, mean, na.rm = TRUE)
df$Lower <- apply(RO, 2, quantile, probs = 0.025, na.rm = TRUE)
df$Upper <- apply(RO, 2, quantile, probs = 0.975, na.rm = TRUE)

# select ages!
plotAges <- c(2, 4, 6, 8, 10, 12)

# plot
df %>%
  mutate(Year = Year + 2007,
         Age  = factor(Age)) %>%
  filter(Age %in% plotAges) %>% 
  ggplot(aes(x = Year, y = Mean, group = Age, colour = Age)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = Age), alpha = 0.2, colour = NA) +
  geom_line(linewidth = 0.8) +
  scale_x_continuous(breaks = c(2008, 2010, 2012, 2014, 2016, 2018, 2020, 2022)) +
  labs(x = "Year", y = "Reproductive output", colour = "Age", fill = "Age") +
  theme_bw()

# ggsave("figures/results12ageCs/routput.jpeg", width = 18.0, height = 10.0, units = c("cm"), dpi = 600)


## Covariate effects on survival -----------------------------------------------

inv_logit <- function(x) 1 / (1 + exp(-x))

keep_ages <- c(1, 10, 11, 12, 13)

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
  mutate(Age = factor(Age)) %>%
  ggplot(aes(x = x, y = Mean, colour = Age)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = Age), alpha = 0.2, colour = NA) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ covariate, scales = "free_x") +
  scale_y_continuous(limits = c(0.2, 1), breaks = c(0.2, 0.4, 0.6, 0.8, 1.0)) +
  labs(x = "Standardized covariate value", y = "Survival",
       colour = "Age", fill = "Age") +
  theme_bw()

# ADD RAW DATA POINTS ON TOP?!

