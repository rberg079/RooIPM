#' Wrangle environmental data
#'
#' @param dens.data # path to xlsx file of population density data to use. As of Apr 2025: "data/WPNP_Methods_Results_January2025.xlsx".
#' @param veg.data # path to xlsx file of vegetation data to use. As of Apr 2025: "data/biomass data April 2009 - Jan 2025_updated Feb2025.xlsx".
#' @param wea.data # path to xlsx file of weather data to use. As of Apr 2025: "data/Prom_Weather_2008-2023_updated Jan2025 RB.xlsx".
#' @param wind.data # path to csv file of wind data to use. As of Apr 2025: "data/POWER_Point_Daily_20080101_20241231_10M.csv".
#'
#' @returns a list containing year, veg, dens, vegE, densE, nNoVeg, & nNoDens
#' @export
#'
#' @examples

wrangleData_env <- function(dens.data, veg.data, wea.data, wind.data){
  
  ## Load libraries
  library(readxl)
  library(tidyverse)
  library(lubridate)
  
  ## Load data
  density <- read_excel(dens.data)
  biomass <- read_excel(veg.data)
  weather <- read_excel(wea.data)
  wind <- read_csv(wind.data, skip = 13)
  
  ## Sort population density data
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
    # approximate sem using CIs (G. Pigeon's code)
    mutate(SE = (DensUCI-DensLCI)/(qnorm(0.975)*2),
           DensSE = SE/sd(Dens, na.rm = T)) %>%
    select(SeasYr, Dens, DensSE)
  
  
  ## Sort vegetation data
  biomass <- biomass %>% 
    rename(Veg = "DW Pal in") %>% 
    select(ID, Day, Month, Year, Veg) %>% 
    mutate(Veg = as.numeric(Veg)*4000, # convert to g/m^2
           Date = ymd(paste(Year, Month, Day, sep = "-")),
           Year = as.numeric(Year),
           Month = month(Date),
           Day = day(Date),
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
  
  # calculate daily vegetation growth per exclosure
  biomass <- biomass %>%
    group_by(ID) %>%
    mutate(Lag = lag(Date),
           Lapse = round(difftime(Date, Lag, units = "days")),
           DailyVeg = as.numeric(Veg) / as.numeric(Lapse)) %>%
    ungroup()
  
  # calculate mean & sd over all exclosures
  biomass <- biomass %>%
    group_by(SeasYr) %>%
    mutate(mDailyVeg = mean(DailyVeg, na.rm = TRUE),
           sdDailyVeg = sd(DailyVeg, na.rm = TRUE)) %>%
    ungroup() %>% 
    filter(!is.na(mDailyVeg)) %>% 
    distinct(Date, mDailyVeg, sdDailyVeg) %>% 
    rename(Veg = "mDailyVeg", VegSE = "sdDailyVeg")
  
  
  ## Sort weather data
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
    select(Date, Year, Month, Day, SeasYr, Rain)
  
  
  ## Sort wind data
  wind <- wind %>%
    rename(Year = "YEAR", Month = "MO", Day = "DY", Max = "T2M_MAX", Min = "T2M_MIN",
           Wind = "WS10M", Gusts = "WS10M_MAX") %>%  # wind is in m/s at 10 m
    mutate(Date = ymd(paste(Year, Month, Day)),
           Wind = -0.863+0.427*Wind,                 # calculate wind at 0.4 m
           Gusts = -0.863+0.427*Gusts) %>%           # calculate wind at 0.4 m
    select(-PRECTOTCORR)
  
  
  ## Join all environmental data &
  ## calculate Nixon-Smith chill index
  env <- weather %>% 
    left_join(biomass) %>% 
    left_join(density) %>%
    fill(Veg, .direction = "up") %>%
    fill(VegSE, .direction = "up") %>% 
    select(Date, Year, Month, Day, SeasYr,
           Dens, DensSE, Veg, VegSE, Rain) %>% 
    mutate(Veg = ifelse(Date < "2009-04-22", NA, Veg),
           VegSE = ifelse(Date < "2009-04-22", NA, VegSE)) %>% 
    left_join(wind)
  
  # Nixon-Smith Chill Index (BOM working paper, 1972)
  # C = (11.7 + 3.1(wind^0.5))(40 - T) + 481 + 418(1 - e^-0.04*R)
  env <- env %>% 
    mutate(Chill = ((11.7 + 3.1*(sqrt(Gusts)))*(40 - Min)) +
             418 + (418*(1 - exp(-0.04*Rain))),
           Warn.18 = ifelse(Chill >= 1000, 1, 0)) %>%  # 0.82 percentile
           # Warn.10 = ifelse(Chill >= 1062, 1, 0),    # 0.90 percentile
           # Warn.05 = ifelse(Chill >= 1124, 1, 0),    # 0.95 percentile
    group_by(Year, Month) %>% 
    mutate(Warns.18 = sum(Warn.18, na.rm = T)) %>% 
           # Warns.05 = sum(Warn.05, na.rm = T),
           # Warns.10 = sum(Warn.10, na.rm = T),
    ungroup() %>% 
    select(Date, Year, Month, Day,
           Dens, DensSE, Veg, VegSE, Rain, Max, Min,
           Wind, Gusts, Chill, Warn.18, Warns.18)
  
  
  ## Summarise by year,
  ## where year X spans Sept 1 X to Aug 31 X+1
  env <- env %>% 
    filter(Date > "2008-03-31") %>%  # 6 mos before 1st census
    mutate(Year = ifelse(Month < 10, Year-1, Year)) %>%
    distinct(Date, Year, Month, Day, Dens, DensSE, Veg, VegSE)
  
  # propagate uncertainty in density & vegetation data
  dens <- env %>% 
    filter(!is.na(Dens)) %>% 
    distinct(Year, Dens, DensSE) %>% 
    mutate(DensSE = DensSE^2) %>% 
    group_by(Year) %>% 
    mutate(Dens = mean(Dens),
           DensSE = sqrt(sum(DensSE)) / 5) %>% 
    ungroup() %>% 
    distinct(Year, Dens, DensSE)
  
  veg <- env %>% 
    filter(!is.na(Veg)) %>% 
    select(Date, Year, Veg, VegSE) %>% 
    group_by(Year) %>% 
    mutate(Veg = sum(Veg),
           VegSE = sqrt(sum(VegSE^2)),
           across(c(Veg, VegSE), ~replace(., Year == 2008, NA)),
           across(c(Veg, VegSE), ~replace(., Year == 2024, NA))) %>% 
    ungroup() %>% 
    distinct(Year, Veg, VegSE)
  
  # join & calculate vegetation per capita
  env <- veg %>% 
    left_join(dens) %>% 
    mutate(VegRoo = Veg / Dens,
           VegRooSE = abs(VegRoo)*sqrt((VegSE / Veg)^2+(DensSE / Dens)^2))
  
  
  ## Return scaled data
  year <- seq(from = 1, to = 17, by = 1)
  
  veg  <- as.numeric(scale(env$Veg))
  dens <- as.numeric(scale(env$Dens))
  
  vegE  <- as.numeric(ifelse(is.na(env$VegSE), 2, env$VegSE/sd(env$Veg, na.rm = T)))
  densE <- as.numeric(ifelse(is.na(env$DensSE), 2, env$DensSE/sd(env$Dens, na.rm = T)))
  
  nNoVeg  <- sum(is.na(veg))
  nNoDens <- sum(is.na(dens))
  
  return(list(year = year,
              veg = veg,
              dens = dens,
              vegE = vegE,
              densE = densE,
              nNoVeg = nNoVeg,
              nNoDens = nNoDens))
  
}


# test <- wrangleData_env(dens.data = "data/WPNP_Methods_Results_January2025.xlsx",
#                         veg.data  = "data/biomass data April 2009 - Jan 2025_updated Feb2025.xlsx",
#                         wea.data  = "data/Prom_Weather_2008-2023_updated Jan2025 RB.xlsx",
#                         wind.data = "data/POWER_Point_Daily_20080101_20241231_10M.csv")
# test

