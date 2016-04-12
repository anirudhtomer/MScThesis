model = function(){
  for(i in 0:(nsubjects-1)){
    for(j in 1:nrep){
      dsweight[i*nrep+j]~dnorm(mu[i*nrep+j],errPrecision)
      mu[i*nrep+j]<- betaGender*dsgender[i*nrep+j] +
        betaBy*dsby[i*nrep+j] +
        randomComp[i+1,1] + randomComp[i+1,2]*dstime[i*nrep+j]
    }
    
    randomComp[i+1,1:2]~dmnorm(randmu[S[i+1], 1:2], randPrecision[1:2,1:2,S[i+1]])
    #randomComp[i+1,1:2]~dmnorm(randmu[S[i+1], 1:2], randPrecision[1:2,1:2])
    S[i+1]~dcat(Eta[])
  }
  
  for(k in 1:ncomponents){
    randmu[k,1]~dnorm(betaMu, betaTau)
    randmu[k,2]~dnorm(betaMu, betaTau)
    
    randSigma[1,1,k]<-1/precision1[k]
    randSigma[1,2,k]<-rho[k] / sqrt(precision1[k] * precision2[k])
    randSigma[2,1,k]<-rho[k] / sqrt(precision1[k] * precision2[k])
    randSigma[2,2,k]<-1/precision2[k]
    rho[k]~dunif(0,1)
    randPrecision[1:2,1:2,k]<-inverse(randSigma[1:2,1:2,k])
    precision1[k]~dgamma(0.0001,0.0001)
    precision2[k]~dgamma(0.0001,0.0001)
    
#     randPrecision[1:2,1:2,k]~dwish(wishartParm[,],3)
#     randSigma[1:2,1:2,k]<-inverse(randPrecision[1:2,1:2,k])
  }
  
  # randPrecision[1:2,1:2]~dwish(wishartParm[,],4)
  # randSigma[1:2,1:2]<-inverse(randPrecision[1:2,1:2])
  
  errPrecision~dgamma(varGammaParm,varGammaParm)
  
  betaGender~dnorm(betaMu,betaTau)
  betaBy~dnorm(betaMu,betaTau)
  
  Eta~ddirch(dirichParm[])
}

singleModel = function(){
  for(i in 0:(nsubjects-1)){
    for(j in 1:nrep){
      dsweight[i*nrep+j]~dnorm(mu[i*nrep+j],errPrecision)
      mu[i*nrep+j]<- betaGender*dsgender[i*nrep+j] +
        betaBy*dsby[i*nrep+j] +
        randomComp[i+1,1] + randomComp[i+1,2]*dstime[i*nrep+j]
    }
    
    randomComp[i+1,1:2]~dmnorm(randmu[S[i+1], 1:2], randPrecision[1:2,1:2,S[i+1]])
    S[i+1]~dcat(Eta[])
  }
  
  randmu[1,1]~dnorm(betaMu, betaTau)
  randmu[1,2]~dnorm(betaMu, betaTau)
  
  randSigma[1,1,1]<-1/precision1[1]
  randSigma[1,2,1]<-rho[1] / sqrt(precision1[1] * precision2[1])
  randSigma[2,1,1]<-rho[1] / sqrt(precision1[1] * precision2[1])
  randSigma[2,2,1]<-1/precision2[1]
  randPrecision[1:2,1:2,1]<-inverse(randSigma[1:2,1:2,1])
  precision1[1]~dgamma(0.0001,0.0001)
  precision2[1]~dgamma(0.0001,0.0001)
  rho[1]~dunif(-1,1)
  
#   randPrecision[1:2,1:2,1]~dwish(wishartParm[,],3)
#   randSigma[1:2,1:2,1]<-inverse(randPrecision[1:2,1:2,1])
  
  errPrecision~dgamma(varGammaParm,varGammaParm)
  
  betaGender~dnorm(betaMu,betaTau)
  betaBy~dnorm(betaMu,betaTau)
  
  Eta[1]~dbeta(100,1)
}

