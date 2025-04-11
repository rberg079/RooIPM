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
  
  
  ## Simulate vital rates ------------------------------------------------------
  
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
  
  
  ## Simulate initial population sizes -----------------------------------------
  
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
  
  
  ## Assemble myinits list -----------------------------------------------------
  
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

