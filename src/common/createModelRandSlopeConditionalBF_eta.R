model_eta = function(){
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
  
  Eta~ddirch(dirichParm[])
}

singleModel_eta = function(){
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
  
  Eta[1]~dbeta(100,1)
}

