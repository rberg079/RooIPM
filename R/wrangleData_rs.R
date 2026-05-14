#' Wrangle reproductive success data
#'
#' @param rs.data character string. Path to xlsx file of reproductive success data to use. As of Apr 2025: "data/RSmainRB_Mar25.xlsx".
#' @param obs.data character string. Path to xlsx file of observation data to use. As of Aug 2025: "data/PromObs_2008-2023.xlsx".
#' @param prime integer vector. Age range which should be considered prime age for reproductive success. prime = c(5:11) by default.
#' @param ageClasses integer. Number of age classes to be considered in the reproductive success model. ageClasses = 20 by default.
#' @param known.age logical. If TRUE, females of unknown age are filtered out. known.age = FALSE by default.
#' @param cum.surv logical. If TRUE, survival is calculated as cumulative survival up to the current step. cum.surv = TRUE by default.
#'
#' @returns a list containing a series of potential response variables, individual and population covariates related to reproductive success.
#' @export
#'
#' @examples

wrangleData_rs <- function(rs.data, obs.data, prime = c(5:11), ageClasses = 20,
                           known.age = FALSE, cum.surv = TRUE){

  # # for testing purposes
  # rs.data = "data/RSmainRB_May26.xlsx"
  # ageClasses = 12
  # known.age = TRUE
  # cum.surv = FALSE
  
  
  ## Set up --------------------------------------------------------------------
  
  # load libraries
  library(readxl)
  suppressPackageStartupMessages(library(tidyverse))
  
  # load data
  rs <- suppressWarnings(read_excel(rs.data))
  
  # clean up
  rs <- suppressWarnings(rs %>% 
    select(ID, Age, Year, Capture, Exclude, Weight, Leg, Teeth, 
           Repro, Parturition, PYid, SurvLPY, SurvWN, SurvNov1, SurvNov2,
           PYLastObs, Dead, HRDead) %>% 
    mutate(Age = as.numeric(Age)) %>% 
    filter(Age <= 20 | is.na(Age), # exclude 21 & 22 year-olds
           Exclude == 0, # exclude new females caught for their young
           # remove first born of "twins", often dropped at March capture
           PYid != 308 & PYid != 340 & PYid != 672 & PYid != 885 & PYid != 900 &
           PYid != 891 & PYid != 912 & PYid != 1023 & PYid != 1106 | is.na(PYid)) %>% 
    # subtract 1 to Year & Age so reproductive attempts are associated with the year
    # in which they began rather than with the cohort year the young is born into
    # so this is all aligned with survival & environmental data
    mutate(Year = Year-1, Age = Age-1,
           # reformat PYLastObs
           PYLastObs = case_when(
             is.na(PYLastObs) ~ NA_Date_,
             TRUE ~ as.Date(paste0("01-", PYLastObs),
                            format = "%d-%m-%Y") %m+% months(1) - days(1)),
           # deduce SurvSep
           SurvSep1 = ifelse(SurvNov1 == 1, 1, NA),
           SurvSep1 = case_when(
             SurvNov1 == 2 ~ NA,
             is.na(SurvSep1) & PYLastObs > as.Date(paste0(Year+1, "-09-01")) ~ 1,
             is.na(SurvSep1) & PYLastObs < as.Date(paste0(Year+1, "-09-01")) ~ 0,
             is.na(SurvSep1) & is.na(PYLastObs) & SurvLPY == 0 ~ 0,
             TRUE ~ SurvSep1),
           SurvSep2 = ifelse(SurvNov2 == 1, 1, NA),
           SurvSep2 = case_when(
             SurvNov2 == 2 ~ NA,
             is.na(SurvSep2) & PYLastObs > as.Date(paste0(Year+2, "-09-01")) ~ 1,
             is.na(SurvSep2) & PYLastObs < as.Date(paste0(Year+2, "-09-01")) ~ 0,
             is.na(SurvSep2) & is.na(PYLastObs) & SurvWN == 0 ~ 0,
             TRUE ~ SurvSep2)))
  
  # limit to females of known age or not
  if(known.age){
    rs <- rs[!is.na(rs$Age),]
  }
  
  
  ## Calculate individual covariates -------------------------------------------

  # age class & birthdate
  rs <- suppressWarnings(rs %>% 
    mutate(Age = as.numeric(Age),
           AgeC = case_when(Age == 0 ~ 1,             # yaf
                            between(Age, 1, 2) ~ 2,   # subadult
                            between(Age, 3, 6) ~ 3,   # prime-age
                            between(Age, 7, 9) ~ 4,   # pre-senescent
                            Age > 10 ~ 5, TRUE ~ NA), # senescent
           Parturition = mdy(Parturition)))
  
  
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
  
  
  ## Return (mostly) scaled data -----------------------------------------------
  
  # birth rate data (B)
  rs <- rs %>% filter(!is.na(Repro))
  
  B <- rs$Repro
  nB <- length(B)
  year.B <- as.integer(factor(rs$Year))
  
  # survival of pouch young data (R)
  rs <- rs %>% filter(!is.na(SurvSep1))
  
  R <- rs$SurvSep1
  nR <- length(R)
  id.R <- as.integer(rs$ID)
  id.R <- match(id.R, sort(unique(id.R)))
  nID.R <- length(unique(id.R))
  year.R <- as.integer(factor(rs$Year))
  
  # sort age & age classes
  # age.R is age in Sept before breeding!
  age.R <- as.integer(rs$Age)
  nAge <- max(age.R)
  
  if(ageClasses == 6){
    ageC.R = c(0,1,2,3,4,4,5,5,5,5, rep(6,30))
  }else if(ageClasses == 12){
    ageC.R = c(0,1,2,3,4,5,6,7,8,9,10, rep(11,29))
  }else if(ageClasses == 20){
    ageC.R = c(seq(from = 0, to = 18, by = 1), rep(18,21))
  }
  nAgeC.R <- max(ageC.R)
  
  return(list(B = B,
              nB = nB,
              year.B = year.B,
              
              R = R,
              nR = nR,
              id.R = id.R,
              nID.R = nID.R,
              year.R = year.R,
              
              age.R = age.R,
              nAge = nAge,
              ageC.R = ageC.R,
              nAgeC.R = nAgeC.R
              ))
  
}

