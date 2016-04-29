hierarchicalModel = function(){
  for(i in 0:(nsubjects-1)){
    for(j in 1:nrep){
      dsweight[i*nrep+j]~dnorm(mu[i*nrep+j],errPrecision)
      mu[i*nrep+j]<- betaGender*dsgender[i*nrep+j] +
        betaBy*dsby[i*nrep+j] +
        betaTime*dstime[i*nrep+j] +
        randomIntercept[i+1]
    }
    randomIntercept[i+1]~dnorm(randmu[S[i+1]],randPrecision[S[i+1]])
    S[i+1]~dcat(Eta[])
  }
  
  for(b in 1:ncomponents){
    randmu[b]~dnorm(betaMu, betaTau)
    randPrecision[b]~dgamma(varGammaParm, varGammaParm)
  }
  
  errPrecision~dgamma(varGammaParm,varGammaParm)
  
  betaGender~dnorm(betaMu,betaTau)
  betaBy~dnorm(betaMu,betaTau)
  betaTime~dnorm(betaMu,betaTau)
  
  Eta~ddirch(dirichParm[])
}

marginalModel = function(){
  for(i in 1:nsubjects){
    weight[i,(1:nrep)]~dmnorm(mu[i,(1:nrep)], omega[1:nrep, 1:nrep,S[i]])
    for(j in 1:nrep){
      XBeta[i,j]<-betaGender*dsgender[i]+betaBy*dsby[i]+betaTime*dstime[j]
    }
    mu[i,(1:nrep)]<-XBeta[i,(1:nrep)] + randmu[S[i]]
    S[i]~dcat(Eta[])
  }
  
  #randmu[1:ncomponents]~dmnorm(betaRandMu[1:ncomponents], temp[1:ncomponents, 1:ncomponents])
  #temp~dwish(betaRandMuTau,ncomponents)
  for(b in 1:ncomponents){
    #mue[b]~dnorm(betaMu, betaTau)
    randmu[b]~dnorm(betaMu, betaTau)
    randPrecision[b]~dgamma(varGammaParm, varGammaParm)
    
    correlation[b]~dunif(0,1)
    
    for(c in 1:(nrep*(nrep-1))){
      sigma[covMatNonDiagIndices[c*2-1],covMatNonDiagIndices[c*2], b] = correlation[b]/randPrecision[b]
    }
    
    for(n in 1:nrep){
      sigma[n,n,b]<-1/randPrecision[b]
    }
    omega[1:nrep, 1:nrep, b]<-inverse(sigma[1:nrep, 1:nrep,b])
  }
  #randmu<-sort(mue)
  
  errPrecision~dgamma(varGammaParm,varGammaParm)
  
  betaGender~dnorm(betaMu,betaTau)
  betaBy~dnorm(betaMu,betaTau)
  betaTime~dnorm(betaMu,betaTau)
  
  Eta~ddirch(dirichParm[])
}

singleMarginalModel = function(){
  for(i in 1:nsubjects){
    weight[i,(1:nrep)]~dmnorm(mu[i,(1:nrep)], omega[1:nrep, 1:nrep,S[i]])
    for(j in 1:nrep){
      XBeta[i,j]<-betaGender*dsgender[i]+betaBy*dsby[i]+betaTime*dstime[j]
    }
    mu[i,(1:nrep)]<-XBeta[i,(1:nrep)] + randmu[S[i]]
    S[i]~dcat(Eta[])
  }
  
  randmu[1]~dnorm(betaMu, betaTau)
  randPrecision[1]~dgamma(varGammaParm, varGammaParm)
    
  for(c in 1:(nrep*(nrep-1))){
    sigma[covMatNonDiagIndices[c*2-1],covMatNonDiagIndices[c*2], 1] = 1/randPrecision[1]
  }
    
  for(n in 1:nrep){
    sigma[n,n,1]<-1/randPrecision[1] + 1/errPrecision
  }
  omega[1:nrep, 1:nrep, 1]<-inverse(sigma[1:nrep, 1:nrep,1])
  
  errPrecision~dgamma(varGammaParm,varGammaParm)
  
  betaGender~dnorm(betaMu,betaTau)
  betaBy~dnorm(betaMu,betaTau)
  betaTime~dnorm(betaMu,betaTau)
  
  Eta[1]~dbeta(100,1)
}
