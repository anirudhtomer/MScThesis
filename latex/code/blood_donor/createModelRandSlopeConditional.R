model = function(){
  for(i in 0:(nsubjects-1)){
    for(j in 1:nrep[[i]]){
      dsHb[i*nrep+j]~dnorm(mu[i*nrep+j],errPrecision)
      mu[i*nrep+j]<- betaDonate*dsDonate[i*nrep+j] + betaAge*dsAge[i*nrep+j] + 
        betaSeason*dsSeason[i*nrep+j] + betaTSPD*dsTSPD[i*nrep+j] +
        randomComp[i+1,1] + randomComp[i+1,2]*(donationLastTwoYears[i*nrep+j]*0.1)
    }
    
    randomComp[i+1,1:2]~dmnorm(randmu[S[i+1], 1:2], randPrecision[1:2,1:2,S[i+1]])
    #randomComp[i+1,1:2]~dmnorm(randmu[S[i+1], 1:2], randPrecision[1:2,1:2])
    S[i+1]~dcat(Eta[])
  }
  
  for(k in 1:ncomponents){
    randmu[k,1]~dnorm(betaMu, betaTau)
    randmu[k,2]~dnorm(betaMu, betaTau)
    #randmu[k,1:2]~dmnorm(betaMuMult, inverse(randSigma[1:2,1:2,k]/10))
    
    # randSigma[1,1,k]<-1/precisionIntercept[k]
    # randSigma[1,2,k]<-rho[k] / sqrt(precisionIntercept[k] * precisionSlope[k])
    # randSigma[2,1,k]<-rho[k] / sqrt(precisionIntercept[k] * precisionSlope[k])
    # randSigma[2,2,k]<-1/precisionSlope[k]
    # rho[k]~dunif(0,1)
    # randPrecision[1:2,1:2,k]<-inverse(randSigma[1:2,1:2,k])
    # precisionIntercept[k]~dgamma(gammaShapeRate, gammaShapeRate)
    # precisionSlope[k]~dgamma(gammaShapeRate, gammaShapeRate)
    
    randPrecision[1:2,1:2,k]~dwish(wishartParm,3)
    randSigma[1:2,1:2,k]<-inverse(randPrecision[1:2,1:2,k])
  }
  
  # randPrecision[1:2,1:2]~dwish(wishartParm[,],4)
  # randSigma[1:2,1:2]<-inverse(randPrecision[1:2,1:2])
  
  errPrecision~dgamma(gammaShapeRate,gammaShapeRate)
  
  betaDonate~dnorm(betaMu,betaTau)
  betaAge~dnorm(betaMu,betaTau)
  betaSeason~dnorm(betaMu,betaTau)
  betaTSPD~dnorm(betaMu,betaTau)
  
  Eta~ddirch(dirichParm[])
}

singleModel = function(){
  for(i in 0:(nsubjects-1)){
    for(j in 1:nrep[[i]]){
      dsHb[i*nrep+j]~dnorm(mu[i*nrep+j],errPrecision)
      mu[i*nrep+j]<- betaDonate*dsDonate[i*nrep+j] + betaAge*dsAge[i*nrep+j] + 
        betaSeason*dsSeason[i*nrep+j] + betaTSPD*dsTSPD[i*nrep+j] +
        randomComp[i+1,1] + randomComp[i+1,2]*(donationLastTwoYears[i*nrep+j]*0.1)
    }
    
    randomComp[i+1,1:2]~dmnorm(randmu[S[i+1], 1:2], randPrecision[1:2,1:2,S[i+1]])
    #randomComp[i+1,1:2]~dmnorm(randmu[S[i+1], 1:2], randPrecision[1:2,1:2])
    S[i+1]~dcat(Eta[])
  }
  
  randmu[1,1]~dnorm(betaMu, betaTau)
  randmu[1,2]~dnorm(betaMu, betaTau)
  
  # randSigma[1,1,1]<-1/precisionIntercept[1]
  # randSigma[1,2,1]<-rho[1] / sqrt(precisionIntercept[1] * precisionSlope[1])
  # randSigma[2,1,1]<-rho[1] / sqrt(precisionIntercept[1] * precisionSlope[1])
  # randSigma[2,2,1]<-1/precisionSlope[1]
  # randPrecision[1:2,1:2,1]<-inverse(randSigma[1:2,1:2,1])
  #precisionIntercept[1]~dgamma(gammaShapeRate,gammaShapeRate)
  #precisionSlope[1]~dgamma(gammaShapeRate,gammaShapeRate)
  #rho[1]~dunif(-1,1)
  
  
  randPrecision[1:2,1:2,1]~dwish(wishartParm, 3)
  randSigma[1:2,1:2,1]<-inverse(randPrecision[1:2,1:2,1])
  
  errPrecision~dgamma(gammaShapeRate,gammaShapeRate)
  
  betaDonate~dnorm(betaMu,betaTau)
  betaAge~dnorm(betaMu,betaTau)
  betaSeason~dnorm(betaMu,betaTau)
  betaTSPD~dnorm(betaMu,betaTau)
  
  Eta[1]~dbeta(100,1)
}

