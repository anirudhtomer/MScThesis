model = function(){
  for(i in 0:(nsubjects-1)){
    for(j in 1:nrep){
      dsweight[i*nrep+j]~dnorm(xbeta[i*nrep+j], errPrecision)
      xbeta[i*nrep+j]<- betaGender*dsgender[i*nrep+j] +
        betaBy*dsby[i*nrep+j] +
        betaTime*dstime[i*nrep+j] +
        randomIntercept[i+1]
    }
    randomIntercept[i+1]~dnorm(randmu[S[i+1]],randPrecision)
    S[i+1]~dcat(Eta[])
  }
  
  for(k in 1:ncomponents){
    mue[k]~dnorm(betaMu, betaTau)
  }
  randmu<-sort(mue)
  
  randPrecision~dgamma(varGammaParm, varGammaParm)
  errPrecision~dgamma(varGammaParm,varGammaParm)
  
  betaGender~dnorm(betaMu,betaTau)
  betaBy~dnorm(betaMu,betaTau)
  betaTime~dnorm(betaMu,betaTau)
  
  Eta~ddirch(dirichParm[])
}

singleModel = function(){
  for(i in 0:(nsubjects-1)){
    for(j in 1:nrep){
      dsweight[i*nrep+j]~dnorm(xbeta[i*nrep+j], errPrecision)
      xbeta[i*nrep+j]<- betaGender*dsgender[i*nrep+j] +
        betaBy*dsby[i*nrep+j] +
        betaTime*dstime[i*nrep+j] +
        randomIntercept[i+1]
    }
    randomIntercept[i+1]~dnorm(randmu[S[i+1]],randPrecision)
  }
  
  randmu[1]~dnorm(betaMu, betaTau)
  
  randPrecision~dgamma(varGammaParm, varGammaParm)
  errPrecision~dgamma(varGammaParm,varGammaParm)
  
  betaGender~dnorm(betaMu,betaTau)
  betaBy~dnorm(betaMu,betaTau)
  betaTime~dnorm(betaMu,betaTau)
  
  Eta[1]~dbeta(100,1)
}