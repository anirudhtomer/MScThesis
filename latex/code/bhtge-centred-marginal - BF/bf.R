modeFunc=function(data){
  d = density(data, n=1e4)
  return(d$x[which.max(d$y)])
}

getEta=function(mcmcfit, summaryFuncName="mean", mcmcIterNum = NA){
  summaryFunc = get(summaryFuncName)
  eta=numeric()
  
  for(i in 1:ncomponents){
    if(!is.numeric(mcmcIterNum)){
      eta[i] = summaryFunc(mcmcfit[[1]][,paste("Eta[",i,"]", sep="")])
    }else{
      eta[i] = mcmcfit[[1]][mcmcIterNum, paste("Eta[",i,"]", sep="")]
    }
  }
  
  return(eta)
}

getRandMu=function(mcmcfit,  summaryFuncName="mean",  mcmcIterNum = NA){
  summaryFunc = get(summaryFuncName)
  randmu=numeric()
  
  for(i in 1:ncomponents){
    if(!is.numeric(mcmcIterNum)){
      randmu[i] = summaryFunc(mcmcfit[[1]][,paste("randmu[",i,"]", sep="")])
    }else{
      randmu[i] = mcmcfit[[1]][mcmcIterNum, paste("randmu[",i,"]", sep="")]
    }
  }
  return(randmu)
}

getSigma=function(mcmcfit,  summaryFuncName="mean", mcmcIterNum = NA){
  summaryFunc = get(summaryFuncName)
  sigma = diag(nrep)
  if(!is.numeric(mcmcIterNum)){
    for(i in 1:nrep){
      for(j in 1:nrep){
        sigma[i,j] = summaryFunc(mcmcfit[[1]][, paste("sigma[",i,",",j,"]", sep="")])
      }
    }
  }else{
    for(i in 1:nrep){
      for(j in 1:nrep){
        sigma[i,j] = mcmcfit[[1]][mcmcIterNum, paste("sigma[",i,",",j,"]", sep="")]
      }
    }
  }
  return(sigma)
}

getOmega=function(mcmcfit,  summaryFuncName="mean", mcmcIterNum = NA){
  summaryFunc = get(summaryFuncName)
  omega = diag(nrep)
  
  if(!is.numeric(mcmcIterNum)){
    for(i in 1:nrep){
      for(j in 1:nrep){
        omega[i,j] = summaryFunc(mcmcfit[[1]][, paste("omega[",i,",",j,"]", sep="")])
      }
    }
  }else{
    for(i in 1:nrep){
      for(j in 1:nrep){
        omega[i,j] = mcmcfit[[1]][mcmcIterNum, paste("omega[",i,",",j,"]", sep="")]
      }
    }
  }
  return(omega)
}

getXBeta=function(mcmcfit,  summaryFuncName="mean", mcmcIterNum = NA){
  summaryFunc = get(summaryFuncName)
  XBeta=matrix(data=NA, nrow=nsubjects, ncol = nrep)
  
  if(!is.numeric(mcmcIterNum)){
    for(i in 1:nsubjects){
      for(j in 1:nrep){
        XBeta[i,j] = summaryFunc(mcmcfit[[1]][,paste("XBeta[",i,",",j,"]", sep="")])
      }
    }
  }else{
    for(i in 1:nsubjects){
      for(j in 1:nrep){
        XBeta[i,j] = mcmcfit[[1]][mcmcIterNum, paste("XBeta[",i,",",j,"]", sep="")]
      }
    }
  }
  return(XBeta)
}


#Step 1: Find theta*

mahalaDist=foreach(i=1:mcmcLen, .combine='rbind') %do%{
 
}
#Step 2: log prior p(theta*)
log_prior=function(eta, randmu, beta, sigma) {
  lgamma(sum(a)) - sum(lgamma(a)) + sum((a-1)*log(w))
  + sum(dnorm(mu, mean = mu0, sd = sqrt(va0), log = TRUE))
  + sum((nu0/2)*log(de0/2) - lgamma(nu0/2) - (nu0/2+1)*log(va) - de0/(2*va))
}

beep(sound=8)
HPDinterval(mcmc(mahalaDist))
quantile(mahalaDist,probs = c(0.025, 0.975))
plot(density(mahalaDist))
