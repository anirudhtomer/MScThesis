model = function(){
  for(i in 1:nsubjects){
    weight[i,(1:nrep)]~dmnorm(mu[i,(1:nrep)], omega)
    for(j in 1:nrep){
      XBeta[i,j]<-betaGender*dsgender[i]+betaBy*dsby[i]+betaTime*dstime[j]
    }
    mu[i,(1:nrep)]<-XBeta[i,(1:nrep)] + randmu[S[i]]
    S[i]~dcat(Eta[])
  }
  
  for(b in 1:ncomponents){
    #mue[b]~dnorm(betaMu, betaTau)
    randmu[b]~dnorm(betaMu, betaTau)
  }
  #randmu<-sort(mue)
  
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
}
