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

getPDandDIC=function(meanPostDeviance, Dthetabar){
  pd = meanPostDeviance - Dthetabar
  dic = meanPostDeviance + pd
  
  list("pd"=pd,"dic"=dic)
}

obsDIC=function(mcmcfit, obsDevParName, summaryFuncName){
  meanPostDeviance = mean(mcmcfit[[1]][, obsDevParName])

  eta = getEta(mcmcfit, summaryFuncName)
  randMu = getRandMu(mcmcfit, summaryFuncName)
  sigma = getSigma(mcmcfit, summaryFuncName)
  omega = getOmega(mcmcfit, summaryFuncName)
  xBeta = getXBeta(mcmcfit, summaryFuncName)
  
  Dthetabar = 0
  obsLikelihood = numeric()
  for(i in 1:nsubjects){
    for(j in 1:ncomponents){
      resid = as.matrix(dswide_y[i,]-xBeta[i,]-randMu[j])
      obsLikelihood[j]<-eta[j]*exp(-0.5 * (t(resid)%*%omega%*%resid))
    }
    
    Dthetabar = Dthetabar + log(det(sigma)) + nrep*log(2*pi) -2*log(sum(obsLikelihood))
  }
  
  getPDandDIC(meanPostDeviance, Dthetabar)
}

calculateObsDIC1 = function(mcmcfit, obsDevParName){
  obsDIC(mcmcfit, obsDevParName, "mean")
}

calculateObsDIC2 = function(mcmcfit, obsDevParName){
  obsDIC(mcmcfit, obsDevParName, "modeFunc")
}

calculateObsDIC3 = function(mcmcfit, obsDevParName, dic3name){
  meanPostDeviance = mean(mcmcfit[[1]][, obsDevParName])
  
  logfhaty  = 0
  for(i in 1:nsubjects){
    logfhaty = logfhaty + log(mean(mcmcfit[[1]][, paste(dic3name,"[",i,"]", sep="")]))
  }
  
  getPDandDIC(meanPostDeviance, -2*logfhaty)
}

calculateDIC7 = function(mcmcfit, obsDevParName){
  meanPostDeviance = mean(mcmcfit[[1]][, obsDevParName])
  
  eta=numeric()
  randmu=numeric()
  for(i in 1:ncomponents){
    eta[i] = mean(mcmcfit[[1]][,paste("eta[",i,"]", sep="")])
    randmu[i] = mean(mcmcfit[[1]][,paste("randmu[",i,"]", sep="")])
  }
  
  precision=mean(mcmcfit[[1]][,"precision"])
  
  Dthetabar = 0
  
  for(i in 1:nobs){
    clusterNum = round(mean(mcmcfit[[1]][,paste("S[",i,"]", sep="")]))
    Dthetabar = Dthetabar + log(sqrt(precision/(2*pi))*exp(-0.5 * ((sample[i]-randmu[clusterNum])^2) * precision))
  }
  
  getPDandDIC(meanPostDeviance, Dthetabar * (-2))
}

calculateCompleteDataMeanPostDeviance=function(mcmcfit){
  mcmcLen = ((niter-nburnin)/nthin)
  totalDev=foreach(i=1:mcmcLen, .combine = c) %dopar%{
    dev = 0
    omega = getOmega(mcmcfit, mcmcIterNum = i)
    eta = getEta(mcmcfit, mcmcIterNum = i)
    for(j in 1:nsubjects){
      clusterNum = round(mcmcfit[[1]][i, paste("S[",j,"]", sep="")])
      
      dev = dev - 2*log(eta[clusterNum])
      dev = dev + log(2*pi) - log(precision)
      dev = dev + precision * ((sample[j]-mcmcfit[[1]][i,paste("mu[",clusterNum,"]",sep="")])^2)
    }
    dev
  }
  sum(totalDev)/mcmcLen
}