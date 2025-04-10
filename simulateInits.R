# 10 April 2025
# Function to simulate initial values for IPM

## Set up ----------------------------------------------------------------------

# load libraries
library(tidyverse)

# load data
# source("wrangleData_rs.R")
# test <- wrangleData_rs(rs.data = "data/RSmainRB_Mar25.xlsx",
#                        obs.data = "data/PromObs_2008-2019.xlsx")

rs <- read_excel("data/RSmainRB_Mar25.xlsx")

# wrangle some RS data
rs <- rs %>% 
  filter(Year > 2007,
         Exclude == 0,
         Age > 2 | is.na(Age),
         between(Year, 2008, 2023)) %>%
  select(ID, Year, Repro, SurvLPY, SurvWN, SurvNov1, SurvNov2, PYLastObs) %>% 
  mutate(SurvLPY = as.numeric(SurvLPY),
         SurvWN = as.numeric(SurvWN),
         # to set SurvLPY & WN at 0 when Repro is 0 instead of NA
         # to represent the probability of producing an LPY/weaning a young
         # as opposed to probability of an existing jellybean surviving to LPY/WN
         SurvLPY = ifelse(Repro == 0, 0, SurvLPY),
         SurvWN = ifelse(Repro == 0, 0, SurvWN),
         # same for SurvNov1 & 2
         SurvNov1 = ifelse(Repro == 0, 0, SurvNov1),
         SurvNov2 = ifelse(Repro == 0, 0, SurvNov2),
         # 2s are specific cases of NA, set to NA
         Repro = ifelse(!is.na(Repro) & Repro == 2, NA, Repro),
         SurvLPY = ifelse(!is.na(SurvLPY) & SurvLPY == 2, NA, SurvLPY),
         SurvWN = ifelse(!is.na(SurvWN) & SurvWN == 2, NA, SurvWN),
         SurvNov1 = ifelse(!is.na(SurvNov1) & SurvNov1 == 2, NA, SurvNov1),
         SurvNov2 = ifelse(!is.na(SurvNov2) & SurvNov2 == 2, NA, SurvNov2),
         # figure out SurvSep from SurvNov & PYLastObs
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


## Simulate vital rates --------------------------------------------------------

# breeding rate
b.exp <- mean(rs$Repro, na.rm = T) # spans 0.58-0.92
b <- runif(1, b.exp*0.75, ifelse(b.exp*1.25 < 1, b.exp*1.25, 1))

# survival of PYs
# (to 1st Sept 1, when they become YAFs)
sPY.exp <- mean(rs$SurvSep1, na.rm = T) # spans 0.32-0.82
s.PY <- runif(1, sPY.exp*0.5, ifelse(sPY.exp*2 < 1, sPY.exp*2, 1))

# survival of YAFs
# (to 2nd Sept 1, when they become SA1s)
# (this would be the 1st age class considered in our published CJS model)
sYAF.exp <- mean(rs$SurvSep2, na.rm = T) # spans 0.32-0.82
s.YAF <- runif(1, sYAF.exp*0.25, ifelse(sYAF.exp*4 < 1, sYAF.exp*4, 1))
# this one maybe should not be based on data, should just be runif(1,0,1)...
# ALTHOUGH as of now, data is prob of weaning a young,
# not of survival from Sep1 to Sep2!!!
# START HERE TOMORROW!




hist(runif(1000, sYAF.exp*0.25, ifelse(sYAF.exp*4 < 1, sYAF.exp*4, 1)))

rs %>% 
  group_by(Year) %>% 
  mutate(mSurvSep2 = mean(SurvSep2, na.rm = T)) %>% 
  distinct(Year, mSurvSep2) %>% 
  ggplot(., aes(x = Year, y = mSurvSep2)) +
  geom_point() + ylim(0,1) + theme_bw()




# to simulate initial values for:
YAF
SA
AD
Ntot

b
s.PY
s.YAF
s.SA[1:2]
s.AD[1:nADs]

