library(mvtnorm)
 
getPDandDIC=function(meanPostDeviance, Dthetabar){
  pd = meanPostDeviance - Dthetabar
  dic = meanPostDeviance + pd
  
  list("pd"=pd,"dic"=dic)
}

calculateObsDThetaBar=function(mcmcfit, summaryFuncName){
  summaryFunc = get(summaryFuncName)
  
  etaSummary = getEta(mcmcfit, summaryFuncName)
  sigmaSummary = getSigma(mcmcfit, summaryFuncName)
  randmuSummary = getRandMu(mcmcfit, summaryFuncName)
  betaBy = summaryFunc(mcmcfit[[1]][,"betaBy"])
  betaGender = summaryFunc(mcmcfit[[1]][,"betaGender"])
  
  Dthetabar = 0
  for(i in 1:nsubjects){
    temp = 0
    for(j in 1:ncomponents){
      randomPartMean = randmuSummary[[j]][1] + randmuSummary[[j]][2]*time
      fixedPartMean = betaBy * (as.numeric(dswide[i, "by"])-1) + 
        betaGender * (as.numeric(dswide[i, "gender"])-1)
      temp = temp + etaSummary[j] * dmvnorm(x = dswide_y[i,], 
                                            mean = randomPartMean + fixedPartMean, 
                                            sigma = sigmaSummary[[j]], log = F)
    }
    Dthetabar = Dthetabar + log(temp)
  }
  Dthetabar = -2*Dthetabar
}

calculateObsMeanPostDeviance = function(mcmcfit){
  
  dswide = dswide
  dswide_y = dswide_y
  time = time
  
  totalDtheta=foreach(k=1:nrow(mcmcfit[[1]]), .combine = c, .packages = c('mvtnorm')) %dopar%{
    
    ncomponents = attributes(mcmcfit)$ncomponents
    nsubjects = attributes(mcmcfit)$nsubjects
    
    source("../common/extractFuncRandSlope.R")
    
    eta = getEta(mcmcfit, mcmcIterNum = k)
    sigma = getSigma(mcmcfit, mcmcIterNum = k)
    randmu = getRandMu(mcmcfit, mcmcIterNum = k)
    betaBy = mcmcfit[[1]][k,"betaBy"]
    betaGender = mcmcfit[[1]][k,"betaGender"]
    
    Dtheta = 0
    for(i in 1:nsubjects){
      temp = 0
      for(j in 1:ncomponents){
        randomPartMean = randmu[[j]][1] + randmu[[j]][2]*time
        fixedPartMean = betaBy * (as.numeric(dswide[i, "by"])-1) + 
          betaGender * (as.numeric(dswide[i, "gender"])-1)
        temp = temp + eta[j] * dmvnorm(x = dswide_y[i,], 
                                       mean = randomPartMean + fixedPartMean, 
                                       sigma = sigma[[j]], log = F)
      }
      Dtheta = Dtheta + log(temp)
    }
    -2*Dtheta
  }
  mean(totalDtheta)
}

calculateObsDIC1 = function(mcmcfit){
  getPDandDIC(calculateObsMeanPostDeviance(mcmcfit), calculateObsDThetaBar(mcmcfit, "mean"))
}

calculateObsDIC2 = function(mcmcfit){
  getPDandDIC(calculateObsMeanPostDeviance(mcmcfit), calculateObsDThetaBar(mcmcfit, "modeFunc"))
}

calculateObsDIC3 = function(mcmcfit){
  
  DThetabar_functionalapprox = 0
  
  dswide = dswide
  dswide_y = dswide_y
  time = time
  obsDev_Subjects=foreach(k=1:nrow(mcmcfit[[1]]), .combine = 'cbind', .packages = c('mvtnorm')) %dopar%{
    
    ncomponents = attributes(mcmcfit)$ncomponents
    nsubjects = attributes(mcmcfit)$nsubjects
    
    source("../common/extractFuncRandSlope.R")
    
    eta = getEta(mcmcfit, mcmcIterNum = k)
    sigma = getSigma(mcmcfit, mcmcIterNum = k)
    randmu = getRandMu(mcmcfit, mcmcIterNum = k)
    betaBy = mcmcfit[[1]][k,"betaBy"]
    betaGender = mcmcfit[[1]][k,"betaGender"]
    
    Dtheta = rep(0, nsubjects)
    for(i in 1:nsubjects){
      for(j in 1:ncomponents){
        randomPartMean = randmu[[j]][1] + randmu[[j]][2]*time
        fixedPartMean = betaBy * (as.numeric(dswide[i, "by"])-1) + 
          betaGender * (as.numeric(dswide[i, "gender"])-1)
        Dtheta[i] = Dtheta[i] + eta[j] * dmvnorm(x = dswide_y[i,], 
                                           mean = randomPartMean + fixedPartMean, 
                                           sigma = sigma[[j]], log = F)
      }
    }
    Dtheta
  }
  
  DThetabar_functionalapprox = -2*sum(log(apply(X = obsDev_Subjects, MARGIN = 1, FUN = mean)))
  
  getPDandDIC(calculateObsMeanPostDeviance(mcmcfit), DThetabar_functionalapprox)
}

calculateDIC7 = function(mcmcfit, obsDevParName){
  meanPostDeviance = 0
  
  dswide = dswide
  dswide_y = dswide_y
  time = time
  foreach(k=1:nrow(mcmcfit[[1]]), .combine = 'c', .packages = c('mvtnorm')) %dopar%{
    ncomponents = attributes(mcmcfit)$ncomponents
    nsubjects = attributes(mcmcfit)$nsubjects
    
    source("../common/extractFuncRandSlope.R")
    
    eta = getEta(mcmcfit, mcmcIterNum = k)
    sigma = getSigma(mcmcfit, mcmcIterNum = k)
    randmu = getRandMu(mcmcfit, mcmcIterNum = k)
    betaBy = mcmcfit[[1]][k,"betaBy"]
    betaGender = mcmcfit[[1]][k,"betaGender"]
    allocations = getAllocations(mcmcfit, mcmcIterNum = k)
    
    for(i in 1:nsubjects){
      
    }
  }
  
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