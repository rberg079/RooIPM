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
  
  
  ## Set up --------------------------------------------------------------------
  
  # load libraries
  library(readxl)
  library(tidyverse)
  
  # load data
  surv <- read_excel(surv.data, sheet = surv.sheet)
  yafs <- read_excel(yafs.data)
  
  # clean up
  surv <- surv %>% 
    arrange(ID) %>% 
    rename(Dead = "Found dead") %>%
    rename_with(~gsub("(\\d{2})to(\\d{2})", "in20\\2", .), everything()) %>% 
    # split "08to09" columns into "in2008" & "in2009" columns
    # 1s in "in20XX" should be read as "seen in Aug-Nov 20XX"
    mutate(in2008 = ifelse(!is.na(in2009), 1, NA)) %>%
    mutate(across(in2009:in2023, ~ifelse(is.na(.x) & !is.na(
      get(sub("\\d{4}", as.integer(gsub("\\D", "", cur_column())) + 1, cur_column()))), 1, .x))) %>%
    select(ID, Sex, Dead, in2008, matches("^in20\\d{2}$"), -in2025, matches("^Age\\d{2}")) %>%
    filter(Sex == 2, ID != 1180, ID != 1183)  # females!
  
  
  ## Sort YAF survival data from RS file ---------------------------------------
  
  # figure out survival to Sept 1
  yafs <- yafs %>% 
    arrange(PYid, Year) %>% 
    filter(PYsex == 2, Exclude == 0, !is.na(PYid)) %>% # females
    select(Year, PYid, Repro, SurvLPY, SurvWN, SurvNov1, SurvNov2, PYLastObs) %>% 
    rename(ID = "PYid") %>% 
    mutate(SurvLPY = as.numeric(SurvLPY),
           SurvWN = as.numeric(SurvWN),
           PYLastObs = case_when(
             is.na(PYLastObs) ~ NA_Date_,
             TRUE ~ as.Date(paste0("01-", PYLastObs),
                            format = "%d-%m-%Y") %m+% months(1) - days(1)),
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
  
  tmp <- yafs %>% 
    filter(SurvSep1 == 1, !is.na(SurvSep2)) %>% 
    select(Year, ID, SurvSep1, SurvSep2)
  
  newbies <- tmp$ID[tmp$SurvSep2 == 0]
  remove(tmp)
  
  yafs <- yafs %>% 
    filter(SurvSep1 == 1, !is.na(SurvSep2)) %>% 
    select(Year, ID, SurvSep1, SurvSep2) %>% 
    pivot_longer(cols = starts_with("SurvSep"), names_to = "Time", values_to = "yafs_val") %>%
    mutate(Year = ifelse(Time == "SurvSep1", Year, Year + 1)) %>%
    complete(ID, Year = seq(from = 2008, to = 2024, by = 1)) %>% 
    select(Year, ID, yafs_val) %>% 
    pivot_wider(names_from = Year,
                values_from = yafs_val,
                names_prefix = "in")
  
  yafs_new <- yafs %>% filter(ID %in% newbies) %>% select(ID, in2008:in2024)
  yafs_old <- yafs %>% filter(!(ID %in% newbies)) %>% select(ID, in2008:in2024)
  
  
  ## Create obs matrix ---------------------------------------------------------
  
  obs <- surv %>% 
    select(ID, in2008:in2024) %>% 
    mutate_at(vars(in2008:in2024), ~ifelse(.> 1 | is.na(.), 0, .))
  
  # add yafs_old
  year_cols <- 2:18
  
  for (i in 1:nrow(yafs_old)) {
    id <- yafs_old$ID[i]
    obs_row <- which(obs$ID == id)
    
    if (length(obs_row) == 1) {
      for (j in year_cols) {
        if (!is.na(yafs_old[i, j]) && yafs_old[i, j] == 1 && obs[obs_row, j] == 0) {
          obs[obs_row, j] <- 1
        }
      }
    }
  }
  
  # add yafs_new
  obs <- rbind(obs, yafs_new) %>% 
    arrange(ID) %>% 
    select(in2008:in2024) %>% 
    mutate(across(everything(), ~ replace_na(., 0)))
  
  # roos observed outside of the study area,
  # roadkilled or poached treated as unobserved
  
  
  ## Create state matrix -------------------------------------------------------
  
  # convert all surv scores to 0 or 1
  state <- surv %>% 
    mutate(HRDead = as.numeric(apply(select(., starts_with("in")) == 2, 1, any)),
           HRDead = replace_na(HRDead, 0)) %>% 
    select(ID, Dead, HRDead, in2008:in2024) %>% 
    mutate_at(vars(in2008:in2024), ~case_when(
      . == 2 ~ 0,  # roadkills are dead
      . == 3 ~ 1,  # observed emigrants are alive
      . == 4 ~ NA, # unobserved emigrants are unknown
      TRUE ~ .))
  
  id <- state %>% select(ID, Dead, HRDead)
  state <- state %>% select(ID, in2008:in2024)
  
  # add yafs_old
  year_cols <- 2:18
  
  for (i in 1:nrow(yafs_old)) {
    who <- yafs_old$ID[i]
    state_row <- which(state$ID == who)
    
    if (length(state_row) == 1) {
      for (j in year_cols) {
        if (!is.na(yafs_old[i, j]) && yafs_old[i, j] == 1 && is.na(state[state_row, j])) {
          state[state_row, j] <- 1
        }
      }
    }
  }
  
  # NAs become 0s for roos found dead...
  state.found <- rbind(state, yafs_new) %>% 
    left_join(id, by = "ID") %>% 
    mutate(Dead = ifelse(ID %in% newbies, 1, Dead)) %>% 
    select(ID, Dead, HRDead, in2008:in2024) %>% 
    filter(!is.na(Dead)) %>% 
    arrange(ID)
  
  for(i in which(is.na(state.found[,20]))){
    state.found[i,(max(which(state.found[i,] == 1)) +1):ncol(state.found)] <- 0
  }
  
  # ...but remain NAs for vanished roos
  state.vanished <- state %>% 
    left_join(id, by = "ID") %>% 
    mutate(Dead = ifelse(ID %in% newbies, 1, Dead)) %>% 
    select(ID, Dead, HRDead, in2008:in2024) %>% 
    filter(is.na(Dead)) %>% 
    mutate_at(vars(in2008:in2024), ~na_if(., 0))
  
  # join all roos again
  state <- bind_rows(state.found, state.vanished) %>% arrange(ID)
  id <- state %>% select(ID, Dead, HRDead)
  
  # roos that were missed one year but seen the next were alive
  state[, 4:20] <- t(apply(state[, 4:20], 1, function(row) {
    ones <- which(row == 1)
    if (length(ones) >= 2) {
      row[min(ones):max(ones)] <- 1
    }
    return(row)
  }))
  
  # create first & last
  state <- state %>% 
    mutate(HRDead = replace_na(HRDead, 0)) %>%
    rowwise() %>% 
    mutate(first = min(which(c_across(4:20) == 1)),
           last = min(which(c_across(4:20) == 0)),
           last = ifelse(HRDead == 1, last-1, last),
           last = ifelse(last == Inf, ncol(obs), last))
  
  # save first & last in id dataframe
  id <- state %>% select(ID, Dead, HRDead, first, last)
  state <- state %>% select(in2008:in2024)
  
  
  ## Create age matrix ---------------------------------------------------------
  
  age <- surv %>% 
    select(ID, Age08:Age24) %>% 
    mutate_at(vars(Age08:Age24), ~ifelse(. == "A", NA, .)) %>% 
    mutate_all(~as.numeric(.)) %>% 
    {names(.)[-1] <- sub("Age", "in20", names(.)[-1]); .}
  
  # get age from YAF data
  yafs_age <- yafs_new %>% 
    mutate(across(in2008:in2024, ~ ifelse(is.na(.), NA, 1 - .)))
  
  # add yafs_new
  age <- rbind(age, yafs_age) %>%
    arrange(ID) %>%
    select(-ID)
  
  # fill in some NAs
  fill_ages <- function(row) {
    if(all(is.na(row))) return(row) # if all NAs, return as is
    first <- which(!is.na(row))[1]  # first non-NA value
    
    # fill forward from first known age
    row[first:length(row)] <- seq(from = row[first], by = 1, length.out = length(row) - first + 1)
    
    # fill backward if needed
    if(first > 1) {
      row[1:(first - 1)] <- seq(from = row[first] - first + 1, by = 1, length.out = first - 1)
    }
    return(row)
  }
  
  age <- as.data.frame(t(apply(age, 1, fill_ages))) %>%
    mutate_all(~replace(., . < 0, NA))
  
  id$uka <- as.logical(is.na(age[,ncol(age)]))
  
  
  ## Return (mostly) scaled data -----------------------------------------------
  
  # remove inds who were only in the dataset 1 year
  noInfo <- id$first == id$last
  length(which(noInfo))
  
  obs   <- unname(as.matrix(obs[!noInfo,]))
  state <- unname(as.matrix(state[!noInfo,]))
  age   <- unname(as.matrix(age[!noInfo,])+1)
  ageC  <- c(1,2,2,3,3,3,3,4,4,4, rep(5,30))
  
  nID   <- nrow(state)
  nYear <- ncol(state)
  
  nAgeC  <- max(ageC, na.rm = T)
  noAge   <- which(is.na(age[,ncol(age)]))
  nNoAge <- length(noAge)
  
  first <- as.numeric(id$first[!noInfo])
  last  <- as.numeric(id$last[!noInfo])
  uka   <- id$uka[!noInfo]
  id    <- as.numeric(id$ID[!noInfo])
  
  return(list(obs = obs,
              state = state,
              age = age,
              ageC = ageC,
              
              nID.sv = nID,
              nYear = nYear,
              
              nAgeC = nAgeC,
              noAge = noAge,
              nNoAge = nNoAge,
              
              first = first,
              last = last,
              uka = uka,
              id = id,
              
              W = diag(nAgeC),
              DF = nAgeC))
  
}

# test <- wrangleData_surv(surv.data = "data/PromSurvivalOct24.xlsx",
#                          yafs.data = "data/RSmainRB_Mar25.xlsx")
# test

