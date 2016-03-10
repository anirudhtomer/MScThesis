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

mahalaDist=foreach(i=1:mcmcLen, .combine='rbind') %do%{
  eta = getEta(mcmcfit, mcmcIterNum = i)
  randmu = getRandMu(mcmcfit, mcmcIterNum = i)
  sigmaEst = getSigma(mcmcfit, mcmcIterNum = i)
  xBeta = getXBeta(mcmcfit, mcmcIterNum = i)
  
  obsMahal = rep(0, nsubjects)
  for(j in 1:nsubjects){
    alloc = mcmcfit[[1]][i, paste("S[",j, "]", sep="")]
    randomPart = dswide_y[j,]-xBeta[j,]
    for(m in 1:ncomponents){
      if(alloc!=m){
        obsMahal[j] = obsMahal[j] + mahalanobis(randomPart, rep(randmu[m], nrep), 
                                                sigmaEst)
      }
    }
    obsMahal[j] = obsMahal[j]/(ncomponents-1)
  }
  
  #newMahal = numeric(0)
  #newObsCount = 0
  #for(alloc in 1:ncomponents){
  #  newObs=mvrnorm(n=round(eta[alloc]*1000), rep(randmu[alloc], nrep), sigmaEst)
   # for(f in 1:nrow(newObs)){
   #   newObsCount = newObsCount + 1
    #  newMahal[newObsCount] = 0
     # for(g in 1:ncomponents){
      #  if(g!=alloc){
       #   newMahal[newObsCount] = newMahal[newObsCount] + 
        #    mahalanobis(newObs[f,], rep(randmu[g], nrep), sigmaEst)
        #}
    #  }
    #  newMahal[newObsCount] = newMahal[newObsCount]/(ncomponents-1)
    #}
  #}

    return(mean(obsMahal))
}
beep(sound=8)
HPDinterval(mcmc(mahalaDist))
quantile(mahalaDist,probs = c(0.025, 0.975))
plot(density(mahalaDist))
