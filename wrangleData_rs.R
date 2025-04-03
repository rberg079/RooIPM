

wrangleData_env <- function(rs.data, obs.data){
  
  ## Load libraries
  library(readxl)
  library(tidyverse)
  
  ## Load data
  rs <- read_excel(rs.data)
  obs <- read_excel(obs.data)
  
  ## Sort reproductive success data
  rs <- rs %>% 
    select(ID, Cohort, Age, Acert, Year, Capture, Exclude, Treatment, 
           Weight, Leg, Teeth, Repro, Parturition, PYsex, PYid, SurvLPY, SurvWN, 
           PYLastObs, PYFound, Dead, HRDead) %>%                                      # REMOVE WHAT IS UNUSED
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
           Teeth = ifelse(between(CaptDay, 182, 366), Teeth, NA),
           Teeth = case_when(Teeth == 0.2 ~ 0.5, Teeth == 0.3 ~ 0.5,
                             Teeth == 0.8 ~ 1.0, Teeth == 1.2 ~ 1.0,
                             Teeth == 1.8 ~ 2.0, TRUE ~ Teeth),
           Teeth = Teeth*2)
  
  ## Calculate individual covariates
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
  
  # body condition & mass gain
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
    mutate(PEff = lag(Eff)) %>%
    ungroup()
  
  ## Sort observation data
  
  ## CONTINUE HERE!!
  
}




## Observation data ------------------------------------------------------------

# sort date & time
obs <- obs %>%
  select(Date, Year, Month, Day, Time, ID, X, Y) %>% 
  mutate(ttime = format(as.POSIXct(Time), format = "%H:%M")) %>% 
  select(-Time) %>% 
  rename(Time = ttime) %>% 
  mutate(X = as.numeric(X),
         Y = as.numeric(Y))

obs <- obs %>%
  filter(X < 40000, X > 32000,              # remove typos in X
         !is.na(ID), !is.na(X), !is.na(X),  # remove NAs in ID, X & Y
         Month >= 7)                        # limit to main field season

# limit to IDs seen at least 10x/year
obs <- obs %>%
  group_by(Year, ID) %>%
  mutate(DaysObs = n_distinct(Date)) %>%
  ungroup() %>% 
  arrange(ID, Year) %>%
  filter(DaysObs >= 10)

# calculate xMed
obs <- obs %>%
  group_by(ID, Year) %>%
  mutate(xMed = median(X, na.rm = T)) %>%
  ungroup()

ggplot(obs, aes(x = X)) + geom_density() + theme_bw()
ggplot(obs, aes(x = xMed)) + geom_density() + theme_bw()

obs <- obs %>%
  distinct(ID, Year, DaysObs, xMed)

df <- left_join(df, obs)


## Environmental data ----------------------------------------------------------

# individual environment!
# from 1 month pre-birth to 7 months post-birth
calculate_Veg <- function(parturition, veg_data) {
  if(is.na(parturition)) return(NA)
  
  start <- parturition %m-% months(1) # 1 month before birth
  end <- parturition %m+% months(7) # 7 months after birth
  window <- veg_data %>% filter(Date >= start & Date <= end)
  
  if(any(is.na(window$Veg))) return(NA)
  
  sum(window$Veg, na.rm = T)
}

calculate_Dens <- function(parturition, dens_data) {
  if(is.na(parturition)) return(NA)
  
  start <- parturition %m-% months(1) # 1 month before birth
  end <- parturition %m+% months(7) # 7 months after birth
  window <- dens_data %>% filter(Date >= start & Date <= end)
  
  if(all(is.na(window$Dens))) return(NA)
  
  mean(window$Dens, na.rm = T)
}

calculate_Win <- function(parturition, win_data) {
  if(is.na(parturition)) return(NA)
  
  start <- parturition %m-% months(1) # 1 month before birth
  end <- parturition %m+% months(7) # 7 months after birth
  window <- win_data %>% filter(Date >= start & Date <= end)
  
  if(all(is.na(window$Warn.18))) return(NA)
  
  sum(window$Warn.18, na.rm = T)
}

df <- df %>%
  rowwise() %>%
  mutate(indVeg7 = calculate_Veg(Parturition, env),
         indDens7 = calculate_Dens(Parturition, env),
         indWin7 = calculate_Win(Parturition, env)) %>%
  ungroup()

# from 8 to 21 months post-birth
calculate_Veg <- function(parturition, veg_data) {
  if(is.na(parturition)) return(NA)
  
  start <- parturition %m+% months(8)
  end <- parturition %m+% months(21)
  window <- veg_data %>% filter(Date >= start & Date <= end)
  
  if(any(is.na(window$Veg))) return(NA)
  
  sum(window$Veg, na.rm = T)
}

calculate_Dens <- function(parturition, dens_data) {
  if(is.na(parturition)) return(NA)
  
  start <- parturition %m+% months(8)
  end <- parturition %m+% months(21)
  window <- dens_data %>% filter(Date >= start & Date <= end)
  
  if(all(is.na(window$Dens))) return(NA)
  
  mean(window$Dens, na.rm = T)
}

calculate_Win <- function(parturition, win_data) {
  if(is.na(parturition)) return(NA)
  
  start <- parturition %m+% months(8)
  end <- parturition %m+% months(21)
  window <- win_data %>% filter(Date >= start & Date <= end)
  
  if(all(is.na(window$Warn.18))) return(NA)
  
  sum(window$Warn.18, na.rm = T)
}

df <- df %>%
  rowwise() %>%
  mutate(indVeg21 = calculate_Veg(Parturition, env),
         indDens21 = calculate_Dens(Parturition, env),
         indWin21 = calculate_Win(Parturition, env)) %>%
  ungroup()


## Population data -------------------------------------------------------------

# mean leg length
# mean body mass
# mean mass gain
# mean condition
df <- df %>%
  group_by(Year) %>%
  mutate(mLeg = mean(Leg, na.rm = T),
         mMass = mean(Mass, na.rm = T),
         mMGain = mean(mGain, na.rm = T),
         mCond = mean(Cond, na.rm = T)) %>%
  ungroup()

# proportion in prime-age
# ratio of young weaned to monitored females
prime <- c(4:9)

df <- df %>% 
  group_by(Year) %>% 
  mutate(nFem = n_distinct(ID),
         nka = sum(!is.na(Age)),                   # number of females of known age
         nSA = sum(SurvWN == 1, na.rm = T),        # number of weaned subadults this cohort
         nPrime = sum(Age %in% prime, na.rm = T),  # number of females of prime age
         PropPrime = nPrime/nka,                   # proportion of prime aged
         Ratio = nSA/nFem) %>%                     # ratio of young weaned
  ungroup()

# previous ratio of young weaned
tmp <- df %>%
  distinct(Year, Ratio) %>%
  mutate(PRatio = lag(Ratio))

df <- left_join(df, tmp)
remove(tmp)

# write csv
# write_csv(df, "RS_Mar25_ind.csv")

