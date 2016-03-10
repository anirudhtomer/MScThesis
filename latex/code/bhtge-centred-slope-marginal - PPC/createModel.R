model = function(){
  for(i in 1:nsubjects){
    weight[i,(1:nrep)]~dmnorm(mu[i,(1:nrep)], omega)
    for(j in 1:nrep){
      #XBeta[i,j]<-betaGender*dsgender[i]+betaBy*dsby[i]+betaTime*dstime[j] + 
       # randmu[S[i], 2]*dstime[j]
      XBeta[i,j]<-betaGender*dsgender[i]+betaBy*dsby[i]+
        randmu[S[i], 2]*dstime[j]
    }
    mu[i,(1:nrep)]<-XBeta[i,(1:nrep)] + randmu[S[i],1]
    S[i]~dcat(Eta[])
  }
  
  for(b in 1:ncomponents){
    #mue[b]~dnorm(betaMu, betaTau)
    randmu[b,1]~dnorm(betaMu, betaTau)
    randmu[b,2]~dnorm(betaMu, betaTau)
  }
  #randmu<-sort(mue)
  
  for(c in 1:(nrep*(nrep-1))){
    sigma[covMatNonDiagIndices[c*2-1], covMatNonDiagIndices[c*2]] = 1/d11Precision +
      (1/d22Precision)*covMatNonDiagIndices[c*2-1]*covMatNonDiagIndices[c*2] +
      d12*(covMatNonDiagIndices[c*2-1] + covMatNonDiagIndices[c*2])
  }
  
  for(n in 1:nrep){
    sigma[n,n]<-(1/d22Precision)*n*n + 2*n*d12 + 1/d11Precision + 1/errPrecision 
  }
  
  d12<-rho12/sqrt(d22Precision * d11Precision)
  
  omega<-inverse(sigma)
  
  rho12~dunif(-1,1)
  d11Precision~dgamma(varGammaParm, varGammaParm)
  d22Precision~dgamma(varGammaParm, varGammaParm)
  errPrecision~dgamma(varGammaParm,varGammaParm)
  
  betaGender~dnorm(betaMu,betaTau)
  betaBy~dnorm(betaMu,betaTau)
  betaTime~dnorm(betaMu,betaTau)
  
  Eta~ddirch(dirichParm[])
}
