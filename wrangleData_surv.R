#' Wrangle survival data
#'
#' @param surv.data path to xlsx file of survival data to use. As of Apr 2025: "data/PromSurvivalOct24.xlsx".
#' @param yafs.data path to xlsx file of RS data to use for survival of YAFs. As of Apr 2025: "data/RSmainRB_Mar25.xlsx".
#' @param surv.sheet sheet to select from survival spreadsheet. surv.sheet = "YEARLY SURV" by default.
#'
#' @returns a list containing the obs, state, & age matrices & other parameters needed for the CJS survival model.
#' @export
#'
#' @examples

wrangleData_surv <- function(surv.data, yafs.data, surv.sheet = "YEARLY SURV"){
  
  ## Load libraries
  library(readxl)
  library(tidyverse)
  
  ## Load data
  surv <- read_excel(surv.data, sheet = surv.sheet)
  yafs <- read_excel(yafs.data)
  
  # clean up
  surv <- surv %>% 
    rename(Dead = "Found dead") %>%
    rename_with(~gsub("(\\d{2})to(\\d{2})", "in20\\2", .), everything()) %>% 
    # split "08to09" columns into "in2008" & "in2009" columns
    # 1s in "in20XX" should be read as "seen in Aug-Nov 20XX"
    mutate(in2008 = ifelse(!is.na(in2009), 1, NA)) %>%
    mutate(across(in2009:in2023, ~ifelse(is.na(.x) & !is.na(
      get(sub("\\d{4}", as.integer(gsub("\\D", "", cur_column())) + 1, cur_column()))), 1, .x))) %>%
    select(ID, Sex, Dead, in2008, matches("^in20\\d{2}$"), -in2025, matches("^Age\\d{2}")) %>%
    filter(Sex == 2, ID != 1180, ID != 1183)  # females!
  
  ## Sort YAF survival data from RS file
  # figure out survival to Oct 1  ...SHOULD BE SEPT 1? TO THINK ABOUT!
  yafs <- yafs %>% 
    filter(PYsex == 2, Exclude == 0, !is.na(PYid)) %>%  # females
    select(Year, PYid, SurvLPY, SurvWN, SurvNov1, SurvNov2, PYLastObs) %>% 
    rename(ID = "PYid") %>% 
    mutate(SurvLPY = as.numeric(SurvLPY),
           SurvWN = as.numeric(SurvWN),
           PYLastObs = case_when(
             is.na(PYLastObs) ~ NA_Date_,
             TRUE ~ as.Date(paste0("01-", PYLastObs),
                            format = "%d-%m-%Y") %m+% months(1) - days(1)),
           SurvOct1 = ifelse(SurvNov1 == 1, 1, NA),
           SurvOct1 = case_when(
             SurvNov1 == 2 ~ NA,
             is.na(SurvOct1) & PYLastObs > as.Date(paste0(Year, "-10-01")) ~ 1,
             is.na(SurvOct1) & PYLastObs < as.Date(paste0(Year, "-10-01")) ~ 0,
             is.na(SurvOct1) & is.na(PYLastObs) & SurvLPY == 0 ~ 0,
             TRUE ~ SurvOct1),
           SurvOct2 = ifelse(SurvNov2 == 1, 1, NA),
           SurvOct2 = case_when(
             SurvNov2 == 2 ~ NA,
             is.na(SurvOct2) & PYLastObs > as.Date(paste0(Year+1, "-10-01")) ~ 1,
             is.na(SurvOct2) & PYLastObs < as.Date(paste0(Year+1, "-10-01")) ~ 0,
             is.na(SurvOct2) & is.na(PYLastObs) & SurvWN == 0 ~ 0,
             TRUE ~ SurvOct2)) %>% 
    filter(SurvOct1 == 1, !is.na(SurvOct2)) %>% 
    select(Year, ID, SurvOct1, SurvOct2) %>% 
    pivot_longer(cols = starts_with("SurvOct"), names_to = "Time", values_to = "yafs_val") %>%
    mutate(Year = ifelse(Time == "SurvOct1", Year, Year + 1)) %>%
    complete(ID, Year = seq(from = 2008, to = 2024, by = 1)) %>% 
    select(Year, ID, yafs_val)
  
  ## Create obs matrix
  # pivot surv data to long format
  obs <- surv %>% 
    select(ID, in2008:in2024) %>% 
    mutate_at(vars(in2008:in2024), ~ifelse(.> 1 | is.na(.), 0, .)) %>% 
    pivot_longer(cols = -ID,
                 names_to = "Year",
                 values_to = "obs_val") %>% 
    mutate(Year = as.integer(sub("in", "", Year)))
  
  # add YAF data
  obs <- full_join(obs, yafs, by = c("ID", "Year"))
  
  obs <- obs %>% 
    mutate(new_val = case_when(
      obs_val == 1 | yafs_val == 1 ~ 1,
      obs_val == 0 & (is.na(yafs_val) | yafs_val == 0) ~ 0,
      TRUE ~ NA_real_
    )) %>% 
    select(ID, Year, new_val) %>% 
    pivot_wider(names_from = Year,
                values_from = new_val,
                names_prefix = "in") %>% 
    select(in2008:in2024)
  
  # roos observed outside of the study area,
  # roadkilled or poached treated as unobserved
  
  ## Create state matrix
  # convert all surv scores to 0 or 1
  state <- surv %>% 
    mutate(HRDead = as.numeric(apply(select(., starts_with("in")) == 2, 1, any)),
           HRDead = replace_na(HRDead, 0)) %>% 
    select(ID, Dead, HRDead, in2008:in2024) %>% 
    mutate_at(vars(in2008:in2024), ~case_when(
      . == 2 ~ 0,   # roadkills are dead
      . == 3 ~ 1,   # observed emigrants are alive
      . == 4 ~ NA,  # unobserved emigrants are unknown
      TRUE ~ .))
  
  # roos missed one year but seen the next were alive
  state <- state %>% 
    mutate(across(in2008:in2023, ~ifelse(. == 0 & !is.na(
      get(sub("\\d{4}", as.integer(gsub("\\D", "", cur_column())) + 1, cur_column()))),
      1, .)))
  
  # NAs become 0s for roos found dead...
  state.found <- state %>% 
    filter(!is.na(Dead))
  
  for(i in which(is.na(state.found[,20]))){
    state.found[i,(max(which(state.found[i,] == 1)) +1):ncol(state.found)] <- 0
  }
  
  # ...but remain NAs for vanished roos
  state.vanished <- state %>% 
    filter(is.na(Dead)) %>% 
    mutate_at(vars(in2008:in2024), ~na_if(., 0))
  
  # join all roos again
  state <- bind_rows(state.found, state.vanished)
  id <- state %>% select(ID, Dead, HRDead)
  
  # add YAF data
  state <- state %>%
    arrange(., ID) %>% 
    select(ID, in2008:in2024) %>% 
    pivot_longer(cols = -ID, 
                 names_to = "Year", 
                 values_to = "state_val") %>%
    mutate(Year = as.integer(sub("in", "", Year)))
  
  state <- full_join(state, yafs, by = c("ID", "Year"))
  
  # create first & last
  state <- state %>% 
    mutate(new_val = case_when(
      state_val == 1 | yafs_val == 1 ~ 1,
      state_val == 0 & (is.na(yafs_val) | yafs_val == 0) ~ 0,
      TRUE ~ NA_real_
    )) %>% 
    select(ID, Year, new_val) %>% 
    pivot_wider(names_from = Year,
                values_from = new_val,
                names_prefix = "in") %>% 
    left_join(id, by = "ID") %>% 
    mutate(HRDead = ifelse(is.na(HRDead), 0, HRDead)) %>% 
    rowwise() %>% 
    mutate(first = min(which(c_across(2:18) == 1)),
           last = min(which(c_across(2:18) == 0)),
           last = ifelse(HRDead == 1, last-1, last),
           last = ifelse(last == Inf, ncol(obs), last))
  
  # save first & last in ID dataframe
  id <- state %>% select(ID, Dead, HRDead, first, last)
  state <- state %>% select(in2008:in2024)
  
  ## Create age matrix
  age <- surv %>% 
    select(ID, Age08:Age24) %>% 
    mutate_at(vars(Age08:Age24), ~ifelse(. == "A", NA, .)) %>% 
    mutate_all(~as.numeric(.)) %>% 
    pivot_longer(cols = -ID,
                 names_to = "Year",
                 values_to = "age_val") %>% 
    mutate(Year = as.integer(sub("Age", "20", Year)))
  
  # get age from YAF data
  yafs <- yafs %>%
    group_by(ID) %>%
    mutate(time = row_number(),
           first = min(time[!is.na(yafs_val)], na.rm = TRUE),
           second = min(time[!is.na(yafs_val) & time > first], na.rm = TRUE),
           yafs_val = case_when(
             time == first ~ 0,
             time == second ~ 1,
             TRUE ~ yafs_val)) %>%
    select(-time, -first, -second) %>%
    ungroup()
  
  # add YAF data
  age <- full_join(age, yafs, by = c("ID", "Year"))
  
  age <- age %>% 
    mutate(new_val = case_when(
      is.na(age_val) & yafs_val == 0 ~ 0,
      is.na(age_val) & yafs_val == 1 ~ 1,
      TRUE ~ age_val
    )) %>% 
    select(ID, Year, new_val) %>% 
    pivot_wider(names_from = Year,
                values_from = new_val,
                names_prefix = "in") %>% 
    select(in2008:in2024)
  
  # fill in some NAs
  fill_ages <- function(row) {
    if(all(is.na(row))) return(row)  # if all NAs, return as is
    first <- which(!is.na(row))[1]   # first non-NA value
    
    # fill forward from first known age
    row[first:length(row)] <- seq(from = row[first], by = 1, length.out = length(row) - first + 1)
    
    # fill backward if needed
    if(first > 1) {
      row[1:(first - 1)] <- seq(from = row[first] - first + 1, by = 1, length.out = first - 1)
    }
    return(row)
  }
  
  age <- as.data.frame(t(apply(age, 1, fill_ages))) %>%
    mutate_all(~ replace(., . < 0, NA))
  
  id$uka <- as.logical(is.na(age[,ncol(age)]))
  
  ## Return (mostly) scaled data
  # remove inds who were only in the dataset 1 year
  noInfo <- id$first == id$last
  # noInfo <- id$last == 1
  
  obs   <- unname(as.matrix(obs[!noInfo,]))
  state <- unname(as.matrix(state[!noInfo,]))
  age   <- unname(as.matrix(age[!noInfo,])+1)
  
  nind   <- nrow(state)
  ntimes <- ncol(state)
  
  ageC   <- c(1,2,2,3,3,3,3,4,4,4, rep(5,50))
  nAge   <- max(ageC, na.rm = T)
  noAge  <- which(is.na(age[,ncol(age)]))
  nNoAge <- length(noAge)
  
  first <- as.numeric(id$first)
  last  <- as.numeric(id$last)
  id <- as.numeric(id$ID[!noInfo])
  
  return(list(obs = obs,
              state = state,
              age = age,
              ageC = ageC,
              nind = nind,
              ntimes = ntimes,
              nAge = nAge,
              noAge = noAge,
              nNoAge = nNoAge,
              first = first,
              last = last,
              W = diag(nAge),
              DF = nAge,
              id = id))
  
}

# test <- wrangleData_surv(surv.data = "data/PromSurvivalOct24.xlsx",
#                          yafs.data = "data/RSmainRB_Mar25.xlsx")
# test

