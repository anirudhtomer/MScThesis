model = function(){
  for(i in 0:(nsubjects-1)){
    for(j in 1:nrep){
      dsweight[i*nrep+j]~dnorm(mu[i*nrep+j],errPrecision)
      mu[i*nrep+j]<- betaGender*dsgender[i*nrep+j] +
        betaBy*dsby[i*nrep+j] +
        randomMu[i+1,1] + randomMu[i+1,2]*dstime[i*nrep+j]
    }
    randomMu[i+1,1:2]~dmnorm(randMuMean[S[i+1], 1:2], randPrecision)
    S[i+1]~dcat(Eta[])
  }
  
  for(k in 1:ncomponents){
    randMuMean[k,1]<-intercept + randMuIntercept[k]
    randMuMean[k,2]<-betaTime + randMuSlope[k]
    randMuIntercept_unordered[k]~dnorm(betaMu, betaTau)
    randMuSlope_unordered[k]~dnorm(betaMu, betaTau)
  }
  
  randMuIntercept<-sort(randMuIntercept_unordered-mean(randMuIntercept_unordered))
  randMuSlope<-sort(randMuSlope_unordered-mean(randMuSlope_unordered))
  
  randPrecision~dwish(wishartParm[,],2)
  randSigma<-inverse(randPrecision)
  errPrecision~dgamma(varGammaParm,varGammaParm)
  
  intercept~dnorm(betaMu,betaTau)
  betaGender~dnorm(betaMu,betaTau)
  betaBy~dnorm(betaMu,betaTau)
  betaTime~dnorm(betaMu,betaTau)
  
  Eta~ddirch(dirichParm[])
}
