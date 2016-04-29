model = function(){
  for(i in 1:nsubjects){
    weight[i,(1:nrep)]~dmnorm(mu[i,(1:nrep)], omega[1:nrep, 1:nrep,S[i]])
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
    
    for(c in 1:(nrep*(nrep-1))){
      sigma[covMatNonDiagIndices[c*2-1], covMatNonDiagIndices[c*2], b] = 1/d11Precision[b] +
        (1/d22Precision[b])*covMatNonDiagIndices[c*2-1]*covMatNonDiagIndices[c*2] +
        d12[b]*(covMatNonDiagIndices[c*2-1] + covMatNonDiagIndices[c*2])
    }
    
    for(n in 1:nrep){
      sigma[n,n,b]<-(1/d22Precision[b])*n*n + 2*n*d12[b] + 1/d11Precision[b] + 1/errPrecision 
    }
    
    d12[b]<-rho12[b]/sqrt(d22Precision[b] * d11Precision[b])
    
    omega[1:nrep, 1:nrep, b]<-inverse(sigma[1:nrep, 1:nrep,b])
    
    rho12[b]~dunif(0,1)
    d11Precision[b]~dgamma(varGammaParm, varGammaParm)
    d22Precision[b]~dgamma(varGammaParm, varGammaParm)
  }
  #randmu<-sort(mue)
  
  errPrecision~dgamma(varGammaParm,varGammaParm)
  
  betaGender~dnorm(betaMu,betaTau)
  betaBy~dnorm(betaMu,betaTau)
  betaTime~dnorm(betaMu,betaTau)
  
  Eta~ddirch(dirichParm[])
}