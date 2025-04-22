# 8 April 2025
# Process model for roo IPM
# 2 censuses per year, Sept & March

# September to March transition (over summer)
# age increases during this transition
for (t in 1:ntimes){
  PY[t] <- b[t] * sum(sepAD[3:nADs,t])
  
  marYAF[t] <- s.sepYAF[t] * sepYAF[t]
  
  marSA[t] <- s.sepSA[1,t] * sepSA[1,t]
  
  marAD[1:2,t] <- 0
  marAD[3,t] <- s.sepSA[2,t] * sepSA[2,t]
  
  for (a in 4:nADs){
    marAD[a,t] <- s.sepAD[a-1,t] * sepAD[a-1,t]
  }
  marNtot[t] <- marYAF[t] + marSA[1:2,t] + sum(marAD[3:nADs,t])
}

# March to September transition (over winter)
# year increases during this transition
for (t in 1:ntimes){
  sepYAF[t+1] <- s.PY[t] * PY[t]
  
  sepSA[1,t+1] <- s.marYAF[t] * marYAF[t]
  sepSA[2,t+1] <- s.marSA[t] * marSA[t]
  
  sepAD[1:2,t+1] <- 0
  
  for (a in 3:nADs){
    sepAD[a,t+1] <- s.marAD[a,t] * marAD[a,t]
  }
  sepNtot[t+1] <- sepYAF[t+1] + sum(sepSA[1:2,t+1]) + sum(sepAD[3:nADs,t+1])
}

# initalize Ns for first Sept (t = 1)
# (will eventually be a parameter)
sepYAF[1] <- 20
sepSA[1:2,1] <- 20
sepAD[1:2,1] <- 0
sepAD[3:nADs,1] <- 10
sepNtot[1] <- sepYAF[1] + sum(sepSA[1:2,1]) + sum(sepAD[3:nADs,1])

