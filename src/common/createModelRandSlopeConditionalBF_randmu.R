model_randmu = function(){
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
  
  for(k in 1:ncomponents){
    randmu_unordered[k,1]~dnorm(betaMu, betaTau)
    randmu_unordered[k,2]~dnorm(betaMu, betaTau)
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

singleModel_randmu = function(){
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
  
  errPrecision~dgamma(gammaShapeRate,gammaShapeRate)
  
  betaGender~dnorm(betaMu,betaTau)
  betaBy~dnorm(betaMu,betaTau)
  betaAge~dnorm(betaMu, betaTau)
  
  Eta[1]~dbeta(100,1)
}

