#' Wrangle environmental data
#'
#' @param dens.data character string. Path to xlsx file of population density data to use. As of Apr 2025: "data/WPNP_Methods_Results_January2025.xlsx".
#' @param veg.data character string. Path to xlsx file of vegetation data to use. As of Apr 2025: "data/biomass data April 2009 - Jan 2025_updated Feb2025.xlsx".
#' @param wea.data character string. Path to xlsx file of weather data to use. As of Apr 2025: "data/Prom_Weather_2008-2023_updated Jan2025 RB.xlsx".
#' @param wind.data character string. Path to csv file of wind data to use. As of Apr 2025: "data/POWER_Point_Daily_20080101_20241231_10M.csv".
#' @param obs.data character string. Path to xlsx file of observation data to use. As of Apr 2025: "data/PromObs_2008-2023.xlsx"
#' @param list.data character string. Path to xlsx file of list of all individuals monitored. As of Apr 2025: "data/PromlistAllNov25.xlsx"
#'
#' @returns a list containing year, ab, dens, veg, win, abE, densE, vegE, nNoDens, nNoVeg, nNoWin, & nNoProp.
#' @export
#'
#' @examples

wrangleData_en <- function(dens.data, veg.data, wea.data, wind.data, obs.data, list.data){
  
  # # for testing purposes
  # # dens.data = "data/abundanceData_Proteus.csv" # OR...
  # dens.data = "data/WPNP_Methods_Results_January2026.xlsx"
  # veg.data  = "data/biomass data April 2009 - July 2025_updated Feb2026.xlsx"
  # wea.data  = "data/Prom_Weather_2008-2023_updated Jan2026 RB.xlsx"
  # wind.data = "data/POWER_Point_Daily_20080101_20260331_10M.csv"
  # obs.data  = "data/PromObs_2008-2024.xlsx"
  # list.data = "data/PromlistAllNov25.xlsx"
  
  
  ## Set up --------------------------------------------------------------------
  
  # load libraries
  library(readxl)
  suppressPackageStartupMessages(library(lubridate))
  suppressPackageStartupMessages(library(tidyverse))
  
  # load data
  # density <- read_csv(dens.data, show_col_types = F) # OR...
  density <- suppressWarnings(read_excel(dens.data))
  biomass <- suppressMessages(read_excel(veg.data))
  weather <- suppressWarnings(read_excel(wea.data))
  wind <- read_csv(wind.data, skip = 13, show_col_types = F)
  obs <- suppressWarnings(read_excel(obs.data))
  list <- suppressMessages(read_excel(list.data))
  
  
  ## Density data --------------------------------------------------------------
  
  # if Dave's data:
  density <- density %>%
    rename(Dens = "Mean density",
           lci = "L 95% CI density",
           uci = "U 95% CI density") %>%
    mutate(SeasYr = paste0(substr(Season, 1, 3), Year),

           # # to approximate SDs from 95% CIs
           # # calculate the scaling factor 'C'
           # C_factor = sqrt(uci / lci),
           # # reverse-engineer the coefficient of variation
           # # CV = sqrt(exp((ln(C) / 1.96)^2) - 1)
           # cv = sqrt(exp((log(C_factor) / 1.96)^2) - 1),
           # 
           # # calculate standard error
           # se_exact = Dens * cv,
           # 
           # # calculate standard deviation
           # DensE = se_exact * sqrt(6))

           # OR set SD = 10% (Heloise's span ~4-7%)
           DensE = 0.1*Dens) %>%
    select(SeasYr, Dens, DensE)
  
  # # if Heloise's:
  # density <- density %>%
  #   rename(Ab = "N",
  #          AbE = "SD",
  #          Dens = "D",
  #          DensE = "SD_D") %>%
  #   mutate(SeasYr = paste0(substr(Season, 1, 3), Year)) %>%
  #   select(SeasYr, Ab, AbE, Dens, DensE) %>%
  #   filter(!is.na(Ab))
  
  
  ## Biomass data --------------------------------------------------------------
  
  biomass <- suppressWarnings(biomass %>% 
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
                           paste(Season, Year, sep = ""))))
  
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
  
  
  ## Weather data --------------------------------------------------------------
  
  weather <- suppressWarnings(weather %>% 
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
    filter(Date > "2007-08-31" & Date < "2026-01-01") %>%
    select(Date, Year, Month, Day, SeasYr, Rain))
  
  
  ## Wind data -----------------------------------------------------------------
  
  wind <- wind %>%
    rename(Year = "YEAR", Month = "MO", Day = "DY", Max = "T2M_MAX", Min = "T2M_MIN",
           Wind = "WS10M", Gusts = "WS10M_MAX") %>%  # wind is in m/s at 10 m
    mutate(Date = ymd(paste(Year, Month, Day)),
           Wind = -0.863+0.427*Wind,                 # calculate wind at 0.4 m
           Gusts = -0.863+0.427*Gusts) %>%           # calculate wind at 0.4 m
    select(-PRECTOTCORR)
  
  
  ## Observation data ----------------------------------------------------------
  
  list <- list %>% 
    rename(ID = "I.D.") %>% 
    select(ID, Sex) %>% 
    filter(Sex == "F" | Sex == "M") %>% 
    mutate(Sex = ifelse(Sex == "F", 1, 0))
    
  obs <- obs %>% 
    select(Date, Year, Month, Day, ID) %>% 
    mutate(Day = yday(Date)) %>% 
    filter(between(Day, 213, 304)) %>% 
    left_join(list, by = "ID", relationship = "many-to-many") %>% 
    group_by(Year) %>% 
    summarise(PropF = mean(Sex == 1, na.rm = T), .groups = "drop") %>% 
    ungroup()
  
  
  ## Join it all ---------------------------------------------------------------
  
  env <- weather %>% 
    left_join(biomass, by = "Date") %>% 
    left_join(density, by = "SeasYr") %>%
    fill(Veg, .direction = "up") %>%
    fill(VegSE, .direction = "up") %>% 
    select(Date, Year, Month, Day, SeasYr,
           Dens, DensE, Veg, VegSE, Rain) %>%
    mutate(Veg = ifelse(Date < "2009-04-22", NA, Veg),
           VegSE = ifelse(Date < "2009-04-22", NA, VegSE)) %>% 
    left_join(wind, by = c("Date", "Year", "Month", "Day"))
  
  # calculate Nixon-Smith Chill Index (BOM working paper, 1972)
  # C = (11.7 + 3.1(wind^0.5))(40 - T) + 481 + 418(1 - e^-0.04*R)
  env <- suppressWarnings(env %>%
    mutate(Chill = ((11.7 + 3.1*(sqrt(Gusts)))*(40 - Min)) + 418 + (418*(1 - exp(-0.04*Rain))),
           Warn.18 = ifelse(Chill >= 1000, 1, 0)) %>%  # 0.82 percentile
           # Warn.10 = ifelse(Chill >= 1062, 1, 0),    # 0.90 percentile
           # Warn.05 = ifelse(Chill >= 1124, 1, 0),    # 0.95 percentile
    group_by(Year, Month) %>% 
    mutate(Warns.18 = sum(Warn.18, na.rm = T)) %>% 
           # Warns.05 = sum(Warn.05, na.rm = T),
           # Warns.10 = sum(Warn.10, na.rm = T),
    ungroup() %>% 
    distinct(Date, Year, Month, Day, SeasYr, Dens, DensE, Veg, VegSE, Warns.18))
             # Rain, Max, Min, Wind, Gusts, Chill, Warn.18, Warns.18
  
  # summarise by year,
  # where year X spans Sept 1 X to Aug 31 X+1
  env <- env %>% 
    filter(Date > "2007-08-31") %>%  # 1st census in Sep 2009
    mutate(Year = ifelse(Month < 9, Year-1, Year))
  
  # propagate uncertainty in density & vegetation data
  dens <- env %>% 
    filter(!is.na(Dens)) %>% 
    filter(grepl("Spr", SeasYr)) %>% 
    distinct(SeasYr, Dens, DensE) %>% 
    mutate(Year = as.integer(str_extract(SeasYr, "\\d{4}"))) %>% 
    select(Year, Dens, DensE) %>% 
    complete(Year = 2007:2025)
  
  veg <- env %>% 
    filter(!is.na(Veg)) %>% 
    select(Date, Year, Veg, VegSE) %>% 
    group_by(Year) %>% 
    mutate(Veg = sum(Veg),
           VegSE = sqrt(sum(VegSE^2)),
           across(c(Veg, VegSE), ~replace(., Year == 2008, NA))) %>% 
    ungroup() %>% 
    distinct(Year, Veg, VegSE) %>% 
    complete(Year = 2007:2024) %>% 
    filter(Year < 2025)
  
  win <- env %>% 
    distinct(Year, Month, Warns.18) %>% 
    group_by(Year) %>% 
    mutate(Win = sum(Warns.18),
           Win = ifelse(Year > 2024, NA, Win)) %>% 
    ungroup() %>% 
    distinct(Year, Win) %>% 
    filter(Year < 2025)
  
  # join & calculate vegetation per capita
  env <- tibble(Year = 2008:2025) %>% 
    left_join(dens, by = "Year") %>% 
    left_join(veg,  by = "Year") %>% 
    left_join(win,  by = "Year") %>% 
    mutate(VegRoo = Veg / Dens,
           VegRooSE = abs(VegRoo) * sqrt((VegSE / Veg)^2 + (DensE / Dens)^2)) %>% 
    left_join(obs, by = "Year")
  
  
  ## Return clean data ---------------------------------------------------------
  
  # centre and scale data
  sc <- function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)
  
  year <- seq(from = 1, to = 17, by = 1)
  
  dens  <- as.numeric(env$Dens)          # length nYear
  veg   <- as.numeric(sc(env$Veg[1:17])) # length nYear-1
  win   <- as.numeric(sc(env$Win[1:17]))  # length nYear-1
  propF <- as.numeric(env$PropF)        # length nYear
  
  densE <- as.numeric(ifelse(is.na(env$DensE), 1, env$DensE))
  vegE  <- as.numeric(ifelse(is.na(env$VegSE[1:17]), 1, env$VegSE[1:17]/sd(env$Veg[1:17], na.rm = T))) # to scale uncertainty too
  
  densM  <- mean(dens, na.rm = T)
  densSD <- sd(dens, na.rm = T)
  
  noDens <- which(is.na(dens))
  noVeg  <- which(is.na(veg))
  noWin  <- which(is.na(win))
  noProp <- which(is.na(propF))
  
  nNoDens <- length(noDens)
  nNoVeg  <- length(noVeg)
  nNoWin  <- length(noWin)
  nNoProp <- length(noProp)
  
  area = rep(76.2, 18)
  
  return(list(year = year,
              area = area,
              dens = dens,
              veg = veg,
              win = win,
              propF = propF,
              densM = densM,
              densSD = densSD,
              densE = densE,
              vegE = vegE,
              noDens = noDens,
              noVeg = noVeg,
              noWin = noWin,
              noProp = noProp,
              nNoDens = nNoDens,
              nNoVeg = nNoVeg,
              nNoWin = nNoWin,
              nNoProp = nNoProp))
  
}

