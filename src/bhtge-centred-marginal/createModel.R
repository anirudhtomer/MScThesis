model = function(){
  for(i in 1:nsubjects){
    weight[i,(1:nrep)]~dmnorm(mu[i,(1:nrep)], inverse(sigma[,]))
    for(j in 1:nrep){
      mu[i,j]<-randmu[S[i]]+betaGender*dsgender[i]+
        betaBy*dsby[i]+betaTime*dstime[j]
    }
    S[i]~dcat(Eta[])
  }
  
  for(k in 1:ncomponents){
    mue[k]~dnorm(betaMu, betaTau)
  }
  randmu<-sort(mue)
  
  for(m in 1:(nrep*(nrep-1))){
    sigma[covMatNonDiagIndices[m*2-1],covMatNonDiagIndices[m*2]] = 1/randPrecision
  }
  
  for(n in 1:nrep){
    sigma[n,n]<-1/randPrecision + 1/errPrecision
  }
  
  randPrecision~dgamma(varGammaParm, varGammaParm)
  errPrecision~dgamma(varGammaParm,varGammaParm)
  
  betaGender~dnorm(betaMu,betaTau)
  betaBy~dnorm(betaMu,betaTau)
  betaTime~dnorm(betaMu,betaTau)
  
  Eta~ddirch(dirichParm[])
}
