#' Wrangle reproductive success data
#'
#' @param rs.data path to xlsx file of reproductive success data to use. As of Apr 2025: "data/RSmainRB_Mar25.xlsx".
#' @param obs.data path to xlsx file of observation data to use. As of Apr 2025: "data/PromObs_2008-2019.xlsx".
#' @param prime age range which should be considered prime age. prime = c(4:9) by default.
#'
#' @returns a list containing a series of individual and population covariates related to reproductive success.
#' @export
#'
#' @examples

wrangleData_rs <- function(rs.data, obs.data, prime = c(4:9)){
  
  
  ## Set up --------------------------------------------------------------------
  
  # load libraries
  library(readxl)
  library(tidyverse)
  
  # load data
  rs <- read_excel(rs.data)
  obs <- read_excel(obs.data)
  
  # clean up
  rs <- rs %>% 
    select(ID, Age, Year, Capture, Exclude, Weight, Leg, Teeth, 
           Repro, Parturition, PYid, SurvLPY, SurvWN, 
           PYLastObs, Dead, HRDead) %>% 
    filter(Exclude == 0,  # exclude new females caught for their young
           # remove first born of "twins", often dropped at March capture
           PYid != 308 & PYid != 340 & PYid != 672 & PYid != 885 & PYid != 900 &
           PYid != 891 & PYid != 912 & PYid != 1023 & PYid != 1106 | is.na(PYid)) %>% 
    rename(Mass = Weight) %>% 
    mutate(Mass = as.numeric(Mass),
           # limit to morphometric data collected during the main field season
           Capture = mdy(Capture),
           CaptDay = yday(Capture),
           Mass = ifelse(between(CaptDay, 182, 366), Mass, NA),
           Leg = ifelse(between(CaptDay, 182, 366), Leg, NA),
           # clean up teeth score so it is an integer between 1 & 6
           Teeth = as.numeric(ifelse(between(CaptDay, 182, 366), Teeth, NA)),
           Teeth = case_when(Teeth == 0.2 ~ 0.5, Teeth == 0.3 ~ 0.5,
                             Teeth == 0.8 ~ 1.0, Teeth == 1.2 ~ 1.0,
                             Teeth == 1.8 ~ 2.0, TRUE ~ Teeth),
           Teeth = Teeth*2)
  
  
  ## Calculate individual covariates -------------------------------------------

  # age class & birthdate
  rs <- rs %>% 
    mutate(Age = as.numeric(Age),
           AgeC = case_when(Age == 0 ~ 1,              # yaf
                            between(Age, 1, 2) ~ 2,    # subadult
                            between(Age, 3, 6) ~ 3,    # prime-age
                            between(Age, 7, 9) ~ 4,    # pre-senescent
                            Age > 10 ~ 5, TRUE ~ NA),  # senescent
           Parturition = mdy(Parturition),
           CohortStart = as.Date(paste(Year-1, "08", "01", sep = "-")),
           CohortDay = as.numeric(difftime(Parturition, CohortStart, units = "days")) + 1,
           PYLastObs = case_when(is.na(PYLastObs) ~ NA_Date_,
                                 TRUE ~ as.Date(paste0("01-", PYLastObs),
                                                format = "%d-%m-%Y") %m+% months(1) - days(1)))
  
  # condition & mass gain
  tmp <- rs %>%
    select(ID, Year, Mass, Leg) %>%
    filter(!is.na(Mass) & !is.na(Leg))
  
  res <- rstandard(lm(log(Mass) ~ log(Leg), data = tmp))
  tmp <- cbind(tmp, res)
  
  rs <- left_join(rs, tmp) %>%
    rename(Cond = res) %>%
    group_by(ID) %>%
    mutate(PCond = lag(Cond),
           PMass = lag(Mass),
           mGain = Mass - PMass) %>%
    ungroup()
  
  remove(tmp)
  
  # previous reproductive success
  rs <- rs %>%
    mutate(Eff = ifelse(SurvWN == 1, 3, NA),
           Eff = ifelse(SurvLPY == 1 & is.na(Eff), 2, Eff),
           Eff = ifelse(Repro == 2, NA,
                        ifelse(Repro == 1 & is.na(Eff), 1,
                        ifelse(Repro == 0, 0, Eff)))) %>%
    arrange(ID, Year) %>%
    group_by(ID) %>%
    mutate(PRS = lag(Eff)) %>%
    ungroup()
  

  ## Sort observation data -----------------------------------------------------
  obs <- obs %>%
    select(Date, Year, Month, Day, Time, ID, X, Y) %>% 
    mutate(ttime = format(as.POSIXct(Time), format = "%H:%M")) %>% 
    select(-Time) %>% 
    rename(Time = ttime) %>% 
    mutate(X = as.numeric(X),
           Y = as.numeric(Y)) %>% 
    filter(X < 40000, X > 32000,              # remove typos in X
           !is.na(ID), !is.na(X), !is.na(X),  # remove NAs in ID, X & Y
           Month >= 7)                        # limit to main field season
  
  # limit to IDs seen at least 10x/year
  # calculate median X coordinate
  obs <- obs %>%
    group_by(Year, ID) %>%
    mutate(DaysObs = n_distinct(Date)) %>%
    ungroup() %>% 
    filter(DaysObs >= 10) %>% 
    group_by(ID, Year) %>%
    mutate(xMed = median(X, na.rm = T)) %>%
    ungroup() %>% 
    distinct(ID, Year, DaysObs, xMed)
  
  rs <- left_join(rs, obs)
  
  
  ## Calculate population covariates -------------------------------------------
  
  # mean leg length, condition, mass & mass gain
  rs <- rs %>%
    group_by(Year) %>%
    mutate(mLeg = mean(Leg, na.rm = T),
           mCond = mean(Cond, na.rm = T),
           mMass = mean(Mass, na.rm = T),
           mMGain = mean(mGain, na.rm = T)) %>%
    ungroup()
  
  # proportion in prime-age
  # ratio of young weaned to monitored females
  rs <- rs %>% 
    group_by(Year) %>% 
    mutate(nFem = n_distinct(ID),                    # number of monitored females
           nKA = sum(!is.na(Age)),                   # number of females of known age
           nSA = sum(SurvWN == 1, na.rm = T),        # number of weaned subadults this cohort
           nPrime = sum(Age %in% prime, na.rm = T),  # number of females of prime age
           pPrime = nPrime/nKA,                      # proportion of prime aged
           Ratio = nSA/nFem) %>%                     # ratio of young weaned
    ungroup()
  
  # previous ratio of young weaned
  tmp <- rs %>%
    distinct(Year, Ratio) %>%
    mutate(PRatio = lag(Ratio))
  
  rs <- left_join(rs, tmp)
  remove(tmp)
  
  
  ## Return (mostly) scaled data -----------------------------------------------
  id <- as.numeric(rs$ID)
  year <- as.numeric(as.factor(rs$Year))
  
  age <- rs$Age  # unscaled!
  teeth <- rs$Teeth  # unscaled!
  leg <- scale(rs$Leg)
  mass <- scale(rs$Mass)
  cond <- scale(rs$Cond)
  prs <- rs$PRS  # unscaled!
  xmed <- scale(rs$xMed)
  
  mcond <- scale(rs$mCond)
  pprime <- scale(rs$pPrime)
  ratio <- scale(rs$Ratio)
  pratio <- scale(rs$PRatio)
  
  return(list(id = id,
              year = year,
              prime = prime,
              age = age,
              teeth = teeth,
              leg = leg,
              mass = mass,
              cond = cond,
              prs = prs,
              xmed = xmed,
              mcond = mcond,
              pprime = pprime,
              ratio = ratio,
              pratio = pratio))
  
}

# test <- wrangleData_rs(rs.data = "data/RSmainRB_Mar25.xlsx",
#                        obs.data = "data/PromObs_2008-2019.xlsx",
#                        prime = c(3:12))
# test

