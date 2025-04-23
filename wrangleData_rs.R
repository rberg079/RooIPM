#' Wrangle reproductive success data
#'
#' @param rs.data path to xlsx file of reproductive success data to use. As of Apr 2025: "data/RSmainRB_Mar25.xlsx".
#' @param obs.data path to xlsx file of observation data to use. As of Apr 2025: "data/PromObs_2008-2019.xlsx".
#' @param prime age range which should be considered prime age. prime = c(4:9) by default.
#' @param known.age should we limit to females of known-age? known.age = FALSE by default.
#' @param cum.surv are we interested in cumulative survival to the current step, rather than survival from the previous step? cum.surv = TRUE by default.
#' @param surv.sep1 are we interested in survival to first September, such that NAs should be removed from SurvSep1? surv.sep1 = FALSE by default.
#' @param surv.sep2 are we interested in survival to second September, such that NAs should be removed from SurvSep2? surv.sep2 = FALSE by default.
#'
#' @returns a list containing a series of potential response variables, individual and population covariates related to reproductive success.
#' @export
#'
#' @examples

wrangleData_rs <- function(rs.data, obs.data, prime = c(4:9),
                           known.age = FALSE, cum.surv = TRUE,
                           surv.sep1 = FALSE, surv.sep2 = FALSE){
  
  
  ## Set up --------------------------------------------------------------------
  
  # # for testing purposes
  # rs.data = "data/RSmainRB_Mar25.xlsx"
  # obs.data = "data/PromObs_2008-2019.xlsx"
  # known.age = TRUE
  # cum.surv = TRUE
  # surv.sep1 = TRUE
  
  # load libraries
  library(readxl)
  library(tidyverse)
  
  # load data
  rs <- read_excel(rs.data)
  # obs <- read_excel(obs.data)
  
  # clean up
  rs <- rs %>% 
    select(ID, Age, Year, Capture, Exclude, Weight, Leg, Teeth, 
           Repro, Parturition, PYid, SurvLPY, SurvWN, SurvNov1, SurvNov2,
           PYLastObs, Dead, HRDead) %>% 
    mutate(Age = as.numeric(Age)) %>% 
    filter(Age <= 20 | is.na(Age),  # exclude 21 & 22 year-olds
           Exclude == 0,  # exclude new females caught for their young
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
           Teeth = ifelse(between(CaptDay, 182, 366), Teeth, NA),
           # clean up teeth score so it is an integer between 1 & 6
           Teeth = as.numeric(Teeth),
           Teeth = case_when(Teeth == 0.2 ~ 0.5, Teeth == 0.3 ~ 0.5,
                             Teeth == 0.8 ~ 1.0, Teeth == 1.2 ~ 1.0,
                             Teeth == 1.8 ~ 2.0, TRUE ~ Teeth),
           # remove a few suspicious ones that are likely typos
           Teeth = case_when(Year == 2008 & ID == 2 ~ NA,
                             Year == 2008 & ID == 21 ~ NA,
                             Year == 2008 & ID == 76 ~ NA,
                             TRUE ~ Teeth),
           Teeth = Teeth*2,
           # figure out SurvSep from SurvNov & PYLastObs
           SurvLPY = as.numeric(SurvLPY),
           SurvWN = as.numeric(SurvWN),
           # reformat PYLastObs
           PYLastObs = case_when(
             is.na(PYLastObs) ~ NA_Date_,
             TRUE ~ as.Date(paste0("01-", PYLastObs),
                            format = "%d-%m-%Y") %m+% months(1) - days(1)),
           # deduce SurvSep
           SurvSep1 = ifelse(SurvNov1 == 1, 1, NA),
           SurvSep1 = case_when(
             SurvNov1 == 2 ~ NA,
             is.na(SurvSep1) & PYLastObs > as.Date(paste0(Year, "-09-01")) ~ 1,
             is.na(SurvSep1) & PYLastObs < as.Date(paste0(Year, "-09-01")) ~ 0,
             is.na(SurvSep1) & is.na(PYLastObs) & SurvLPY == 0 ~ 0,
             TRUE ~ SurvSep1),
           SurvSep2 = ifelse(SurvNov2 == 1, 1, NA),
           SurvSep2 = case_when(
             SurvNov2 == 2 ~ NA,
             is.na(SurvSep2) & PYLastObs > as.Date(paste0(Year+1, "-09-01")) ~ 1,
             is.na(SurvSep2) & PYLastObs < as.Date(paste0(Year+1, "-09-01")) ~ 0,
             is.na(SurvSep2) & is.na(PYLastObs) & SurvWN == 0 ~ 0,
             TRUE ~ SurvSep2))
  
  # limit to females of known age or not
  if(known.age){
    rs <- rs[!is.na(rs$Age),]
  }
  
  
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
  
  # # condition & mass gain
  # tmp <- rs %>%
  #   select(ID, Year, Mass, Leg) %>%
  #   filter(!is.na(Mass) & !is.na(Leg))
  # 
  # res <- rstandard(lm(log(Mass) ~ log(Leg), data = tmp))
  # tmp <- cbind(tmp, res)
  # 
  # rs <- left_join(rs, tmp) %>%
  #   rename(Cond = res) %>%
  #   group_by(ID) %>%
  #   mutate(PCond = lag(Cond),
  #          PMass = lag(Mass),
  #          mGain = Mass - PMass) %>%
  #   ungroup()
  # 
  # remove(tmp)
  # 
  # # previous reproductive success
  # rs <- rs %>%
  #   mutate(Eff = ifelse(SurvWN == 1, 3, NA),
  #          Eff = ifelse(SurvLPY == 1 & is.na(Eff), 2, Eff),
  #          Eff = ifelse(Repro == 2, NA,
  #                       ifelse(Repro == 1 & is.na(Eff), 1,
  #                       ifelse(Repro == 0, 0, Eff)))) %>%
  #   arrange(ID, Year) %>%
  #   group_by(ID) %>%
  #   mutate(PRS = lag(Eff)) %>%
  #   ungroup()
  # 
  # 
  # ## Sort observation data -----------------------------------------------------
  # obs <- obs %>%
  #   select(Date, Year, Month, Day, Time, ID, X, Y) %>% 
  #   mutate(ttime = format(as.POSIXct(Time), format = "%H:%M")) %>% 
  #   select(-Time) %>% 
  #   rename(Time = ttime) %>% 
  #   mutate(X = as.numeric(X),
  #          Y = as.numeric(Y)) %>% 
  #   filter(X < 40000, X > 32000,              # remove typos in X
  #          !is.na(ID), !is.na(X), !is.na(X),  # remove NAs in ID, X & Y
  #          Month >= 7)                        # limit to main field season
  # 
  # # limit to IDs seen at least 10x/year
  # # calculate median X coordinate
  # obs <- obs %>%
  #   group_by(Year, ID) %>%
  #   mutate(DaysObs = n_distinct(Date)) %>%
  #   ungroup() %>% 
  #   filter(DaysObs >= 10) %>% 
  #   group_by(ID, Year) %>%
  #   mutate(xMed = median(X, na.rm = T)) %>%
  #   ungroup() %>% 
  #   distinct(ID, Year, DaysObs, xMed)
  # 
  # rs <- left_join(rs, obs)
  # 
  # 
  # ## Calculate population covariates -------------------------------------------
  # 
  # # mean leg length, condition, mass & mass gain
  # rs <- rs %>%
  #   group_by(Year) %>%
  #   mutate(mLeg = mean(Leg, na.rm = T),
  #          mCond = mean(Cond, na.rm = T),
  #          mMass = mean(Mass, na.rm = T),
  #          mMGain = mean(mGain, na.rm = T)) %>%
  #   ungroup()
  # 
  # # proportion in prime-age
  # # ratio of young weaned to monitored females
  # rs <- rs %>% 
  #   group_by(Year) %>% 
  #   mutate(nFem = n_distinct(ID),                    # number of monitored females
  #          nKA = sum(!is.na(Age)),                   # number of females of known age
  #          nSA = sum(SurvWN == 1, na.rm = T),        # number of weaned subadults this cohort
  #          nPrime = sum(Age %in% prime, na.rm = T),  # number of females of prime age
  #          pPrime = nPrime/nKA,                      # proportion of prime aged
  #          Ratio = nSA/nFem) %>%                     # ratio of young weaned
  #   ungroup()
  # 
  # # previous ratio of young weaned
  # tmp <- rs %>%
  #   distinct(Year, Ratio) %>%
  #   mutate(PRatio = lag(Ratio))
  # 
  # rs <- left_join(rs, tmp)
  # remove(tmp)
  
  
  ## Go through toggles --------------------------------------------------------
  
  # if cum.surv is TRUE
  # we want to represent cumulative survival
  # meaning survival up to the current stage (LPY, WN, Sep1, Sep2...)
  # rather than survival from the previous stage to the current stage
  if(cum.surv){
    rs <- rs %>% 
      # ...SurvLPY & WN should thus be 0 when Repro is 0
      mutate(SurvLPY = ifelse(Repro == 0, 0, SurvLPY),
             SurvWN = ifelse(Repro == 0, 0, SurvWN))
  }else{
    rs <- rs %>% 
      # SurvLPY & WN are by default both ONLY conditional on Repro
      # if cum.surv is F, then SurvWN should be NA when SurvLPY is 0
      mutate(SurvWN = ifelse(SurvLPY == 0 | is.na(SurvLPY), NA, SurvWN))
  }
  
  # in all cases
  # 2 are a particular case of NA
  rs <- rs %>% 
    mutate(Repro = ifelse(!is.na(Repro) & Repro == 2, NA, Repro),
           SurvLPY = ifelse(!is.na(SurvLPY) & SurvLPY == 2, NA, SurvLPY),
           SurvWN = ifelse(!is.na(SurvWN) & SurvWN == 2, NA, SurvWN))
  
  # if surv.sep1 is TRUE
  # we are interested in survival to first September
  # & therefore want to remove NAs from SurvSep1
  if(surv.sep1){
    rs <- rs %>% filter(!is.na(SurvSep1))
  }
  
  # if surv.sep2 is TRUE
  # we are interested in survival to second September
  # & therefore want to remove NAs from SurvSep2
  if(surv.sep2){
    rs <- rs %>% filter(!is.na(SurvSep2))
  }
  
  
  ## Return (mostly) scaled data -----------------------------------------------
  
  id   <- as.integer(rs$ID)
  id   <- match(id, sort(unique(id)))
  year <- as.integer(factor(rs$Year))
  
  age  <- as.integer(rs$Age)
  ageC <- c(1,2,2,3,3,3,3,4,4,4, rep(5,30))
  
  N      <- length(id)
  N.id   <- length(unique(id))
  N.year <- length(unique(year))
  N.age  <- length(unique(age))
  N.ageC <- length(unique(ageC))
  
  # teeth <- rs$Teeth       # unscaled!
  # leg <- scale(rs$Leg)
  # mass <- scale(rs$Mass)
  # cond <- scale(rs$Cond)
  # prs <- rs$PRS           # unscaled!
  # xmed <- scale(rs$xMed)
  
  surv7 <- rs$SurvLPY
  surv21 <- rs$SurvWN
  survS1 <- rs$SurvSep1
  survS2 <- rs$SurvSep2
  
  # mcond <- scale(rs$mCond)
  # pprime <- scale(rs$pPrime)
  # ratio <- scale(rs$Ratio)
  # pratio <- scale(rs$PRatio)
  
  return(list(N = N,
              N.id = N.id,
              N.year = N.year,
              N.age = N.age,
              N.ageC = N.ageC,
              id = id,
              year = year,
              age = age,
              ageC = ageC,
              # teeth = teeth,
              # leg = leg,
              # mass = mass,
              # cond = cond,
              # prs = prs,
              # xmed = xmed,
              surv7 = surv7,
              surv21 = surv21,
              survS1 = survS1,
              survS2 = survS2
              # mcond = mcond,
              # pprime = pprime,
              # ratio = ratio,
              # pratio = pratio
              ))
  
}

# test <- wrangleData_rs(rs.data = "data/RSmainRB_Mar25.xlsx",
#                        obs.data = "data/PromObs_2008-2019.xlsx",
#                        prime = c(3:12), known.age = TRUE, cum.surv = TRUE,
#                        surv.sep1 = TRUE, surv.sep2 = FALSE)
# test

