source("../common/createModelRandSlopeConditionalBF_eta.R")
source("../common/createModelRandSlopeConditionalBF_beta.R")
source("../common/createModelRandSlopeConditionalBF_errPrecision.R")
source("../common/createModelRandSlopeConditionalBF_randmu.R")

model = function(){
  for(i in 0:(nsubjects-1)){
    for(j in 1:nrep){
      dsweight[i*nrep+j]~dnorm(mu[i*nrep+j],errPrecision)
      mu[i*nrep+j]<- betaGender*dsgender[i*nrep+j] +
        betaBy*dsby[i*nrep+j] + betaAge*dsage[i*nrep+j]+
        randomComp[i+1,1] + randomComp[i+1,2]*dstime[i*nrep+j]
    }
    
    randomComp[i+1,1:2]~dmnorm(randmu[S[i+1], 1:2], randPrecision[1:2,1:2,S[i+1]])
    #randomComp[i+1,1:2]~dmnorm(randmu[S[i+1], 1:2], randPrecision[1:2,1:2])
    S[i+1]~dcat(Eta[])
  }
  
  for(k in 1:ncomponents){
    randmu_unordered[k,1]~dnorm(betaMu, betaTau)
    randmu_unordered[k,2]~dnorm(betaMu, betaTau)
    # randSigma[1,1,k]<-1/precisionIntercept[k]
    # randSigma[1,2,k]<-rho[k] / sqrt(precisionIntercept[k] * precisionSlope[k])
    # randSigma[2,1,k]<-rho[k] / sqrt(precisionIntercept[k] * precisionSlope[k])
    # randSigma[2,2,k]<-1/precisionSlope[k]
    # rho[k]~dunif(-1,1)
    # randPrecision[1:2,1:2,k]<-inverse(randSigma[1:2,1:2,k])
    # precisionIntercept[k]~dgamma(gammaShapeRate, gammaShapeRate)
    # precisionSlope[k]~dgamma(gammaShapeRate, gammaShapeRate)
    
    randPrecision[1:2,1:2,k]~dwish(wishartPriorScale, wishartPriorDf)
    randSigma[1:2,1:2,k]<-inverse(randPrecision[1:2,1:2,k])
    precisionIntercept[k]<-randPrecision[1,1,k]
    precisionSlope[k]<-randPrecision[2,2,k]
    rho[k]<-randSigma[1,2,k]* sqrt(precisionIntercept[k] * precisionSlope[k])
  }
  orderRandMu<-order(randmu_unordered[1:ncomponents,1])
  for(h in 1:ncomponents){
    randmu[h,1:2]<-randmu_unordered[orderRandMu[h],1:2]
  }
  
  errPrecision~dgamma(gammaShapeRate,gammaShapeRate)
  
  betaGender~dnorm(betaMu,betaTau)
  betaBy~dnorm(betaMu,betaTau)
  betaAge~dnorm(betaMu, betaTau)
  
  Eta~ddirch(dirichParm[])
}

singleModel = function(){
  for(i in 0:(nsubjects-1)){
    for(j in 1:nrep){
      dsweight[i*nrep+j]~dnorm(mu[i*nrep+j],errPrecision)
      mu[i*nrep+j]<- betaGender*dsgender[i*nrep+j] +
        betaBy*dsby[i*nrep+j] + betaAge*dsage[i*nrep+j]+
        randomComp[i+1,1] + randomComp[i+1,2]*dstime[i*nrep+j]
    }
    
    randomComp[i+1,1:2]~dmnorm(randmu[S[i+1], 1:2], randPrecision[1:2,1:2,S[i+1]])
    S[i+1]~dcat(Eta[])
  }
  
  randmu[1,1]~dnorm(betaMu, betaTau)
  randmu[1,2]~dnorm(betaMu, betaTau)
  
  # randSigma[1,1,1]<-1/precisionIntercept[1]
  # randSigma[1,2,1]<-rho[1] / sqrt(precisionIntercept[1] * precisionSlope[1])
  # randSigma[2,1,1]<-rho[1] / sqrt(precisionIntercept[1] * precisionSlope[1])
  # randSigma[2,2,1]<-1/precisionSlope[1]
  # randPrecision[1:2,1:2,1]<-inverse(randSigma[1:2,1:2,1])
  # precisionIntercept[1]~dgamma(gammaShapeRate,gammaShapeRate)
  # precisionSlope[1]~dgamma(gammaShapeRate,gammaShapeRate)
  # rho[1]~dunif(-1,1)
  
  randPrecision[1:2,1:2,1]~dwish(wishartPriorScale, wishartPriorDf)
  randSigma[1:2,1:2,1]<-inverse(randPrecision[1:2,1:2,1])
  precisionIntercept[1]<-randPrecision[1,1,1]
  precisionSlope[1]<-randPrecision[2,2,1]
  rho[1]<-randSigma[1,2,1]* sqrt(precisionIntercept[1] * precisionSlope[1])
  
  errPrecision~dgamma(gammaShapeRate,gammaShapeRate)
  
  betaGender~dnorm(betaMu,betaTau)
  betaBy~dnorm(betaMu,betaTau)
  betaAge~dnorm(betaMu, betaTau)
  
  Eta[1]~dbeta(100,1)
}

