model = function(){
  for(i in 1:nsubjects){
    weight[i,(1:nrep)]~dmnorm(mu[i,(1:nrep)], omega)
    for(j in 1:nrep){
      XBeta[i,j]<-betaGender*dsgender[i]+betaBy*dsby[i]+betaTime*dstime[j]
    }
    
    mu[i,(1:nrep)]<-XBeta[i,(1:nrep)] + randmu[S[i]]
    
    for(a in 1:ncomponents){
      obsLikelihood[i,a]<-Eta[a]*
        exp(-0.5 * t(weight[i,(1:nrep)] - XBeta[i,(1:nrep)] - randmu[a])%*%omega%*%(weight[i,(1:nrep)] - XBeta[i,(1:nrep)] - randmu[a]))
    }
    
    S[i]~dcat(Eta[])
    
    obsLl[i] <- nrep*log(2*pi) + logdet(sigma) -2*log(sum(obsLikelihood[i,]))
    obsL[i]<-exp(obsLl[i]/(-2))
    
    regLl[i] <- nrep*log(2*pi) + logdet(sigma) + 
      t(weight[i,(1:nrep)] - mu[i,(1:nrep)])%*%omega%*%(weight[i,(1:nrep)] - mu[i,(1:nrep)])
  }
  
  for(b in 1:ncomponents){
    mue[b]~dnorm(betaMu, betaTau)
  }
  randmu<-sort(mue)
  
  for(c in 1:(nrep*(nrep-1))){
    sigma[covMatNonDiagIndices[c*2-1],covMatNonDiagIndices[c*2]] = 1/randPrecision
  }
  
  for(n in 1:nrep){
    sigma[n,n]<-1/randPrecision + 1/errPrecision
  }
  
  omega<-inverse(sigma)
  
  randPrecision~dgamma(varGammaParm, varGammaParm)
  errPrecision~dgamma(varGammaParm,varGammaParm)
  
  betaGender~dnorm(betaMu,betaTau)
  betaBy~dnorm(betaMu,betaTau)
  betaTime~dnorm(betaMu,betaTau)
  
  Eta~ddirch(dirichParm[])
  
  regDeviance <- sum(regLl)
  obsDeviance <- sum(obsLl)
}
