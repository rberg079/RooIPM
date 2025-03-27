# 19 March 2025
# Prep environmental data
# Population density, biomass, & weather

library(readxl)
library(tidyverse)
library(lubridate)

# setwd("C:/Users/rberg/OneDrive/Bureau/IPM analyses/Roo IPM")

## Load data -------------------------------------------------------------------

density <- read_excel("WPNP_Methods_Results_January2025.xlsx")
biomass <- read_excel("biomass data April 2009 - Jan 2025_updated Feb2025.xlsx")
weather <- read_excel("Prom_Weather_2008-2023_updated Jan2025 RB.xlsx")
wind <- read_csv("POWER_Point_Daily_20080101_20231231_10M.csv", skip = 13)


## Density ---------------------------------------------------------------------

# Select relevant variables
# Create Year & Season variables
density <- density %>%
  rename(Date = "Month / year",
         Dens = "Mean density",
         DensLCI = "L 95% CI density",
         DensUCI = "U 95% CI density",
         Area = "Area sampled (ha)") %>% 
  mutate(Date = ymd(Date),
         Year = year(Date),
         Month = month(Date),
         Season = case_when(Month <= 2  ~ "Sum",  # Jan & Feb
                            Month <= 5  ~ "Aut",  # Mar to May
                            Month <= 8  ~ "Win",  # Jun to Aug
                            Month <= 11 ~ "Spr",  # Sep to Nov
                            Month == 12 ~ "Sum",  # Dec
                            TRUE ~ NA_character_),
         NextYr = Year +1,
         SeasYr = ifelse(Month == 12,
                         paste(Season, NextYr, sep = ""),
                         paste(Season, Year, sep = ""))) %>% 
  filter(Date > "2008-11-30" & Date < "2025-03-01") %>% 
  select(SeasYr, Dens, DensLCI, DensUCI)


## Biomass ---------------------------------------------------------------------

# Select relevant variables
# Create Date, Year & Season variables
biomass <- biomass %>% 
  rename(Veg = "DW Pal in") %>% 
  select(ID, Day, Month, Year, Veg) %>% 
  mutate(Veg = as.numeric(Veg)*4000, # kg/quadrat to g/m^2
         Date = ymd(paste(Year, Month, Day, sep = "-")),
         Year = as.numeric(Year),
         Month = month(Date),
         Day = as.numeric(Year),
         Season = case_when(Month <= 2  ~ "Sum",  # Jan & Feb
                            Month <= 5  ~ "Aut",  # Mar to May
                            Month <= 8  ~ "Win",  # Jun to Aug
                            Month <= 11 ~ "Spr",  # Sep to Nov
                            Month == 12 ~ "Sum",  # Dec
                            TRUE ~ NA_character_),
         NextYr = as.numeric(Year) +1,
         SeasYr = ifelse(Month == 12,
                         paste(Season, NextYr, sep = ""),
                         paste(Season, Year, sep = "")))

# Calculate daily biomass growth per exclosure
biomass <- biomass %>%
  group_by(ID) %>%
  mutate(Lag = lag(Date),
         Lapse = round(difftime(Date, Lag, units = "days")),
         DailyVeg = as.numeric(Veg) / as.numeric(Lapse)) %>%
  ungroup()

# Calculate mean & sd over all exclosures
biomass <- biomass %>%
  group_by(SeasYr) %>%
  mutate(mDailyVeg = mean(DailyVeg, na.rm = TRUE),
         sdDailyVeg = sd(DailyVeg, na.rm = TRUE)) %>%
  ungroup() %>% 
  filter(!is.na(mDailyVeg)) %>% 
  distinct(Date, mDailyVeg, sdDailyVeg)


## Weather ---------------------------------------------------------------------

# Select relevant variables
# Create Date, Year & Season variables
weather <- weather %>% 
  select(Year, Month, Day, Rain) %>% 
  mutate(Date = ymd(paste(Year, Month, Day, sep = "-")),
         Month = month(Date),
         Season = case_when(Month <= 2  ~ "Sum",  # Jan & Feb
                            Month <= 5  ~ "Aut",  # Mar to May
                            Month <= 8  ~ "Win",  # Jun to Aug
                            Month <= 11 ~ "Spr",  # Sep to Nov
                            Month == 12 ~ "Sum",  # Dec
                            TRUE ~ NA_character_),
         NextYr = as.numeric(Year) +1,
         SeasYr = ifelse(Month == 12,
                         paste(Season, NextYr, sep = ""),
                         paste(Season, Year, sep = ""))) %>%
  filter(Date > "2007-07-31" & Date < "2025-03-01") %>%
  select(Date, Year, Month, Day, Season, SeasYr, Rain)


## Wind ------------------------------------------------------------------------

# Join weather, biomass & density datasets
env <- weather %>% 
  left_join(biomass) %>% 
  left_join(density) %>%
  fill(mDailyVeg, .direction = "up") %>%
  fill(sdDailyVeg, .direction = "up")

# Rename some things
# Replace filled Veg values with NAs before July 2009
env <- env %>% 
  rename(Veg = mDailyVeg, VegSD = sdDailyVeg) %>% 
  select(1:6, Dens, DensLCI, DensUCI, Veg, VegSD, Rain) %>% 
  mutate(Veg = ifelse(Date < "2009-07-30", NA, Veg),
         VegSD = ifelse(Date < "2009-07-30", NA, VegSD))

# Prep & join wind data
wind <- wind %>%
  rename(Year = "YEAR", Month = "MO", Day = "DY", Max = "T2M_MAX", Min = "T2M_MIN",
         Wind = "WS10M", Gusts = "WS10M_MAX") %>%  # wind is in m/s at 10 m
  mutate(Date = ymd(paste(Year, Month, Day)),
         Wind = -0.863+0.427*Wind,                 # calculate wind at 0.4 m
         Gusts = -0.863+0.427*Gusts) %>%           # calculate wind at 0.4 m
  select(-PRECTOTCORR)

env <- env %>% left_join(wind)

# Nixon-Smith Chill Index (BOM working paper, 1972)
# C = (11.7 + 3.1(wind^0.5))(40 - T) + 481 + 418(1 - e^-0.04*R)
env <- env %>% 
  mutate(A = (11.7 + 3.1*(sqrt(Gusts)))*(40 - Min),
         B = 418*(1 - exp(-0.04*Rain)),
         C = A + 481 + B) %>% 
  rename(Chill = C)

quantile(env$Chill, probs = c(0.05, 0.10, 0.25, 0.5, 0.75, 0.90, 0.95), na.rm = T)
percentile <- ecdf(env$Chill)
percentile(1000)

# warnings at >1000 means on 0.18 of days
# warnings on 0.10 of days means at >1062
# warnings on 0.05 of days means at >1124

env <- env %>% 
  mutate(Warn.05 = ifelse(Chill >= 1124, 1, 0),      # 0.95 percentile
         Warn.10 = ifelse(Chill >= 1062, 1, 0),      # 0.90 percentile
         Warn.18 = ifelse(Chill >= 1000, 1, 0)) %>%  # 0.82 percentile
  group_by(Year, Month) %>% 
  mutate(Warns.05 = sum(Warn.05, na.rm = T),
         Warns.10 = sum(Warn.10, na.rm = T),
         Warns.18 = sum(Warn.18, na.rm = T)) %>% 
  ungroup()

# Final variable selection
env <- env %>% 
  select(Date, Year, Month, Day,
         Dens, DensLCI, DensUCI, Veg, VegSD,
         Rain, Max, Min, Wind, Gusts, Chill, Warn.18, Warns.18)

# write csv
# write_csv(env, "Env_Mar25.csv")

