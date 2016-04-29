modeFunc=function(data){
  d = density(data, n=1e4)
  return(d$x[which.max(d$y)])
}

obsDIC=function(mcmcfit, obsDevParName, summaryFuncName){
  meanPostDeviance = mean(mcmcfit[[1]][, obsDevParName])
  
  summaryFunc = get(summaryFuncName)
  
  eta=numeric()
  mu=numeric()
  for(i in 1:ncomponents){
    eta[i] = summaryFunc(mcmcfit[[1]][,paste("eta[",i,"]", sep="")])
    mu[i] = summaryFunc(mcmcfit[[1]][,paste("mu[",i,"]", sep="")])
  }
  
  precision=summaryFunc(mcmcfit[[1]][,"precision"])
  
  Dthetabar = 0
  ObsLikelihood = numeric()
  for(i in 1:nobs){
    for(m in 1:ncomponents){
      ObsLikelihood[m]<-eta[m]*sqrt(precision/(2*pi))*exp(-0.5 * ((sample[i]-mu[m])^2) * precision)
    }
    
    Dthetabar = Dthetabar + log(sum(ObsLikelihood))
  }
  
  getPDandDIC(meanPostDeviance, Dthetabar * (-2))
}

getPDandDIC=function(meanPostDeviance, Dthetabar){
  pd = meanPostDeviance - Dthetabar
  dic = meanPostDeviance + pd
  
  list("pd"=pd,"dic"=dic)
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
  for(i in 1:nobs){
    logfhaty = logfhaty + log(mean(mcmcfit[[1]][, paste(dic3name,"[",i,"]", sep="")]))
  }
  
  getPDandDIC(meanPostDeviance, -2*logfhaty)
}

calculateDIC7 = function(mcmcfit, obsDevParName){
  meanPostDeviance = mean(mcmcfit[[1]][, obsDevParName])
  
  eta=numeric()
  mu=numeric()
  for(i in 1:ncomponents){
    eta[i] = mean(mcmcfit[[1]][,paste("eta[",i,"]", sep="")])
    mu[i] = mean(mcmcfit[[1]][,paste("mu[",i,"]", sep="")])
  }
  
  precision=mean(mcmcfit[[1]][,"precision"])
  
  Dthetabar = 0
  
  for(i in 1:nobs){
    clusterNum = round(mean(mcmcfit[[1]][,paste("S[",i,"]", sep="")]))
    Dthetabar = Dthetabar + log(sqrt(precision/(2*pi))*exp(-0.5 * ((sample[i]-mu[clusterNum])^2) * precision))
  }
  
  getPDandDIC(meanPostDeviance, Dthetabar * (-2))
}