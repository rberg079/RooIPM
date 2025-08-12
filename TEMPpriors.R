# 12 August 2025
# IPM Population model priors to try
# assigning different age-class-specific estimates to true ages

# priors

# survival by age
# from estimates by age class
if(nAgeC.S == 6){
  for(t in 1:(nYear-1)){
    sYF[t] <- S[1, t]
    sSA[t] <- S[2, t]
    sAD[1, t] <- 0
    sAD[2, t] <- S[3, t]
    for(a in 3:6) sAD[a, t] <- S[4, t] # prime-aged
    for(a in 7:9) sAD[a, t] <- S[5, t] # pre-senescent
    for(a in 10:nAge) sAD[a, t] <- S[6, t] # senescent
  }
}else if(nAgeC.S == 12){
  for(t in 1:(nYear-1)){
    sYF[t] <- S[1, t]
    sSA[t] <- S[2, t]
    sAD[1, t] <- 0
    for(a in 2:11) sAD[a, t] <- S[a+1, t] # other adults
    for(a in 12:nAge) sAD[a, t] <- S[13, t] # greybeards
  }
}else if(nAgeC.S == 19){
  for(t in 1:(nYear-1)){
    sYF[t] <- S[1, t]
    sSA[t] <- S[2, t]
    sAD[1, t] <- 0
    for(a in 2:19) sAD[a, t] <- S[a+1, t] # adults
  }
}

# nAgeC.R == 6
ageC.R = c(1,2,3,4,4,4,5,5,5,5, 6, 6, 6, 6, 6, 6, 6, 6, 6 )
trueAR = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19)

# nAgeC.R == 12
ageC.R = c(1,2,3,4,5,6,7,8,9,10,11,12,12,12,12,12,12,12,12)
trueAR = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19)

# reproductive success by age
# from estimates by age class
if(nAgeC.R == 6){
  
  for(t in 1:(nYear-1)){
    
  }
  
}else if(nAgeC.R == 12){
  
}else if(nAgeC.R == 19){
  
}


for(t in 1:(nYear-1)){

  for(a in 1:nAgeC.R){
    sPY[a, 1:(nYear-1)] <- Ra[a, 1:(nYear-1)]
  }
  
  if(nAge > nAgeC.R){
    for(a in (nAgeC.R+1):nAge){
      sPY[a, 1:(nYear-1)] <- Ra[nAgeC.R, 1:(nYear-1)]
    }
  }
}