# 8 April 2025
# Start process model for roo IPM

# September to March transition (over summer)
# age increases during this transition
for (t in 1:ntimes){
  
  PY[t] <- pB[t] * sum(sepAD[3:nADs,t])
  
  marYAF[t] <- s_sepYAF[t] * sepYAF[t]
    
  marSA[t] <- s_sepSA[1,t] * sepSA[1,t]
  
  marAD[1:2,t] <- 0
  
  marAD[3,t] <- s_sepSA[2,t] * sepSA[2,t]
  
  for (a in 4:nADs){
    
    marAD[a,t] <- s_sepAD[a-1,t] * sepAD[a-1,t]
    
  }
  
  marNtot[t] <- marSA[t] + sum(marAD[3:nADs,t])

}

# March to September transition (over winter)
# year increases during this transition
for (t in 1:ntimes){
  
  sepYAF[t+1] <- s_PY[t] * PY[t]
  
  sepSA[1,t+1] <- s_marYAF[t] * marYAF[t]
  
  sepSA[2,t+1] <- s_marSA[t] * marSA[t]
  
  sepAD[1:2,t+1] <- 0
  
  for (a in 3:nADs){
    
    sepAD[a,t+1] <- s_marAD[a,t] * marAD[a,t]
    
  }
  
  sepNtot[t+1] <- sepSA[t+1] + sum(sepAD[3:nADs,t+1])
  
}

# initalize Ns for first Sept (t = 1)
# (will eventually be a parameter)
sepYAF[1] <- 20

sepSA[1,1] <- 20

sepSA[2,1] <- 20

sepAD[1:2,1] <- 0

sepAD[3:nADs,1] <- 10

sepNtot[1] <- sepSA[1] + sum(sepAD[3:nADs,1])



