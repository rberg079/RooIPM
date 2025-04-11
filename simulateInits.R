#' Simulate initial values for IPM
#'
#' @param ntimes integer. Number of time steps in the model. ntimes = 17 by default.
#' @param nADs integer. Number of adult age classes, or maximum age. nADs = 22 by default.
#'
#' @returns list containing initial values for b, s.PY, s.YAF, s.SA, s.AD, YAF, SA, AD & Ntot.
#' @export
#'
#' @examples

simulateInits <- function(ntimes = 17, nADs = 22){
  
  ## Set up ----------------------------------------------------------------------
  
  # rs <- read_excel("data/RSmainRB_Mar25.xlsx")
  # 
  # # wrangle some RS data
  # # to base simulated initial values on
  # # (or at least get inspiration from) real data
  # rs <- rs %>% 
  #   filter(Year > 2007,
  #          Exclude == 0,
  #          Age > 2 | is.na(Age),
  #          between(Year, 2008, 2023)) %>%
  #   select(ID, Year, Repro, SurvLPY, SurvWN, SurvNov1, SurvNov2, PYLastObs) %>% 
  #   mutate(SurvLPY = as.numeric(SurvLPY),
  #          SurvWN = as.numeric(SurvWN),
  #          # # to set SurvLPY & WN at 0 when Repro is 0 instead of NA
  #          # # to represent the probability of producing an LPY/weaning a young
  #          # # as opposed to probability of an existing jellybean surviving to LPY/WN
  #          # SurvLPY = ifelse(Repro == 0, 0, SurvLPY),
  #          # SurvWN = ifelse(Repro == 0, 0, SurvWN),
  #          # # same for SurvNov1 & 2
  #          # SurvNov1 = ifelse(Repro == 0, 0, SurvNov1),
  #          # SurvNov2 = ifelse(Repro == 0, 0, SurvNov2),
  #          # 2s are specific cases of NA, set to NA
  #          Repro = ifelse(!is.na(Repro) & Repro == 2, NA, Repro),
  #          SurvLPY = ifelse(!is.na(SurvLPY) & SurvLPY == 2, NA, SurvLPY),
  #          SurvWN = ifelse(!is.na(SurvWN) & SurvWN == 2, NA, SurvWN),
  #          SurvNov1 = ifelse(!is.na(SurvNov1) & SurvNov1 == 2, NA, SurvNov1),
  #          SurvNov2 = ifelse(!is.na(SurvNov2) & SurvNov2 == 2, NA, SurvNov2),
  #          # figure out SurvSep from SurvNov & PYLastObs
  #          PYLastObs = case_when(
  #            is.na(PYLastObs) ~ NA_Date_,
  #            TRUE ~ as.Date(paste0("01-", PYLastObs),
  #                           format = "%d-%m-%Y") %m+% months(1) - days(1)),
  #          SurvSep1 = ifelse(SurvNov1 == 1, 1, NA),
  #          SurvSep1 = case_when(
  #            SurvNov1 == 2 ~ NA,
  #            is.na(SurvSep1) & PYLastObs > as.Date(paste0(Year, "-09-01")) ~ 1,
  #            is.na(SurvSep1) & PYLastObs < as.Date(paste0(Year, "-09-01")) ~ 0,
  #            is.na(SurvSep1) & is.na(PYLastObs) & SurvLPY == 0 ~ 0,
  #            TRUE ~ SurvSep1),
  #          SurvSep2 = ifelse(SurvNov2 == 1, 1, NA),
  #          SurvSep2 = case_when(
  #            SurvNov2 == 2 ~ NA,
  #            is.na(SurvSep2) & PYLastObs > as.Date(paste0(Year+1, "-09-01")) ~ 1,
  #            is.na(SurvSep2) & PYLastObs < as.Date(paste0(Year+1, "-09-01")) ~ 0,
  #            is.na(SurvSep2) & is.na(PYLastObs) & SurvWN == 0 ~ 0,
  #            TRUE ~ SurvSep2))
  
  
  ## Simulate vital rates --------------------------------------------------------
  
  # breeding rate
  # b.exp <- mean(rs$Repro, na.rm = T) # spans 0.58-0.92
  # b <- runif(1, b.exp*0.75, ifelse(b.exp*1.25 < 1, b.exp*1.25, 1))
  b <- runif(ntimes, 0.5, 1)
  
  # survival of PYs
  # to 1st Sept 1 when they become YAFs
  # sPY.exp <- mean(rs$SurvSep1, na.rm = T) # spans 0.24-0.95
  # s.PY <- runif(1, sPY.exp*0.5, ifelse(sPY.exp*2 < 1, sPY.exp*2, 1))
  s.PY <- runif(ntimes, 0.1, 1)
  
  # survival of YAFs
  # to 2nd Sept 1 when they become SA1s
  # 1st age class considered in our published CJS model
  # sYAF.exp <- mean(rs$SurvSep2, na.rm = T) # spans 0.01-0.88
  s.YAF <- runif(ntimes, 0, 1)
  
  # survival of SA1s to SA2 & SA2 to AD3
  # 2nd age class in our published CJS model
  s.SA <- matrix(runif(2 * (ntimes), 0.5, 1), nrow = 2)
  
  # survival of all ADs
  s.AD            <- matrix(0, nrow = nADs, ncol = ntimes)          # 0s for AD1 & AD2
  s.AD[3:6, ]     <- matrix(runif(4 * (ntimes), 0.6, 1), nrow = 4)  # Prime-age
  s.AD[7:9, ]     <- matrix(runif(3 * (ntimes), 0.5, 1), nrow = 3)  # Pre-senescent
  s.AD[10:nADs, ] <- matrix(runif((nADs - 9) * (ntimes), 0.1, 1),   # Senescent
                            nrow = nADs - 9)
  
  
  ## Simulate initial population sizes -------------------------------------------
  
  # Actual numbers in 2008:
  # 5 female YAFs in Sept, 6 SA1s, 5 SA2s, 21 adults
  # Wendy estimated 22.6% of the population was marked
  
  YAF    <- c(5*5, rep(NA, times = ntimes-1))
  SA     <- matrix(NA, nrow = 2, ncol = ntimes)
  SA[,1] <- c(6*5, 5*5)
  
  AD     <- matrix(NA, nrow = nADs, ncol = ntimes)
  AD[,1] <- c(0, 0, rep(2*5, times = 8), rep(1*5, times = 8))
  
  Ntot   <- c(YAF[1] + sum(SA[1:2,1]) + sum(AD[3:nADs,1]),
              rep(NA, times = ntimes-1))
  
  
  ## Assemble myinits list -------------------------------------------------------
  
  return(list(YAF = YAF,
              SA = SA,
              AD = AD,
              Ntot = Ntot,
              b = b,
              s.PY = s.PY,
              s.YAF = s.YAF,
              s.SA = s.SA,
              s.AD = s.AD))
  
}

# test <- simulateInits(ntimes = 40, nADs = 20)
# test

