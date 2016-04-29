model = function(){
  for(i in 0:(nsubjects-1)){
    for(j in 1:nrep){
      dsweight[i*nrep+j]~dnorm(XBeta[i*nrep+j], errPrecision)
      XBeta[i*nrep+j]<- betaGender*dsgender[i*nrep+j] +
        betaBy*dsby[i*nrep+j] +
        betaTime*dstime[i*nrep+j] +
        randomIntercept[i+1]
    }
    randomIntercept[i+1]~dnorm(randmu[S[i+1]], randPrecision[S[i+1]])
    S[i+1]~dcat(Eta[])
  }
  
  for(k in 1:ncomponents){
    randmue[k]~dnorm(betaMu, betaTau)
    randPrecision[k]~dgamma(gammaShape, gammaRate)
  }
  randmu <- sort(randmue)
  
  errPrecision~dgamma(gammaShape, gammaRate)
  
  betaGender~dnorm(betaMu,betaTau)
  betaBy~dnorm(betaMu,betaTau)
  betaTime~dnorm(betaMu,betaTau)
  
  Eta~ddirch(dirichParm[])
}

####### Model for randmu #######
model_randmu = function(){
  for(i in 0:(nsubjects-1)){
    for(j in 1:nrep){
      dsweight[i*nrep+j]~dnorm(XBeta[i*nrep+j], errPrecision)
      XBeta[i*nrep+j]<- betaGender*dsgender[i*nrep+j] +
        betaBy*dsby[i*nrep+j] +
        betaTime*dstime[i*nrep+j] +
        randomIntercept[i+1]
    }
    randomIntercept[i+1]~dnorm(randmu[S[i+1]], randPrecision[S[i+1]])
    S[i+1]~dcat(Eta[])
  }
  
  for(k in 1:ncomponents){
    randmue[k]~dnorm(betaMu, betaTau)
  }
  randmu <- sort(randmue)
  
  errPrecision~dgamma(gammaShape, gammaRate)
  
  betaGender~dnorm(betaMu,betaTau)
  betaBy~dnorm(betaMu,betaTau)
  betaTime~dnorm(betaMu,betaTau)
  
  Eta~ddirch(dirichParm[])
}

model_errPrecision=function(){
  for(i in 0:(nsubjects-1)){
    for(j in 1:nrep){
      dsweight[i*nrep+j]~dnorm(XBeta[i*nrep+j], errPrecision)
      XBeta[i*nrep+j]<- betaGender*dsgender[i*nrep+j] +
        betaBy*dsby[i*nrep+j] +
        betaTime*dstime[i*nrep+j] +
        randomIntercept[i+1]
    }
    randomIntercept[i+1]~dnorm(randmu[S[i+1]], randPrecision[S[i+1]])
    S[i+1]~dcat(Eta[])
  }
  
  errPrecision~dgamma(gammaShape, gammaRate)
  
  betaGender~dnorm(betaMu,betaTau)
  betaBy~dnorm(betaMu,betaTau)
  betaTime~dnorm(betaMu,betaTau)
  
  Eta~ddirch(dirichParm[])
}

model_beta=function(){
  for(i in 0:(nsubjects-1)){
    for(j in 1:nrep){
      dsweight[i*nrep+j]~dnorm(XBeta[i*nrep+j], errPrecision)
      XBeta[i*nrep+j]<- betaGender*dsgender[i*nrep+j] +
        betaBy*dsby[i*nrep+j] +
        betaTime*dstime[i*nrep+j] +
        randomIntercept[i+1]
    }
    randomIntercept[i+1]~dnorm(randmu[S[i+1]], randPrecision[S[i+1]])
    S[i+1]~dcat(Eta[])
  }
  
  betaGender~dnorm(betaMu,betaTau)
  betaBy~dnorm(betaMu,betaTau)
  betaTime~dnorm(betaMu,betaTau)
  
  Eta~ddirch(dirichParm[])
}

model_eta=function(){
  for(i in 0:(nsubjects-1)){
    for(j in 1:nrep){
      dsweight[i*nrep+j]~dnorm(XBeta[i*nrep+j], errPrecision)
      XBeta[i*nrep+j]<- betaGender*dsgender[i*nrep+j] +
        betaBy*dsby[i*nrep+j] +
        betaTime*dstime[i*nrep+j] +
        randomIntercept[i+1]
    }
    randomIntercept[i+1]~dnorm(randmu[S[i+1]], randPrecision[S[i+1]])
    S[i+1]~dcat(Eta[])
  }
  
  Eta~ddirch(dirichParm[])
}