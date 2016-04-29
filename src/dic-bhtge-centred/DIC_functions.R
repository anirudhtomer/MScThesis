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
  sigma = matrix(nrow=nrep, ncol = nrep)
  
  if(!is.numeric(mcmcIterNum)){
    summaryFunc = get(summaryFuncName)
    errPrecision = summaryFunc(mcmcfit[[1]][, "errPrecision"])
    randPrecision = summaryFunc(mcmcfit[[1]][, "randPrecision"])
    
    for(i in 1:nrep){
      for(j in 1:nrep){
        if(i==j){
          sigma[i,j] = 1/errPrecision + 1/randPrecision
        }else{
          sigma[i,j] = 1/randPrecision
        }
      }
    }
  }else{
    errPrecision = mcmcfit[[1]][i, "errPrecision"]
    randPrecision = mcmcfit[[1]][i, "randPrecision"]
    
    for(i in 1:nrep){
      for(j in 1:nrep){
        if(i==j){
          sigma[i,j] = 1/errPrecision + 1/randPrecision
        }else{
          sigma[i,j] = 1/randPrecision
        }
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

getObsLogL=function(nsubjects,nrep, ncomponents,dsweight, betaTime,dstime,betaBy,dsby,
                    betaGender, dsgender, randMu, eta, detSigma, omega){
  logL = 0
  for(j in 1:nsubjects){
    likelihood = 0
    beginIndex = (j-1)*nrep + 1
    endIndex = j*nrep
    for(k in 1:ncomponents){
      Y = dsweight[beginIndex:endIndex]
      XBeta = betaTime*dstime[beginIndex:endIndex] + betaBy*dsby[beginIndex:endIndex] 
      + betaGender*dsgender[beginIndex:endIndex]
      bias = matrix(Y-XBeta-randMu[k])
      
      likelihood = likelihood + eta[k]*exp(-0.5*(t(bias)%*%omega%*%bias))
    }
    logL = logL - 0.5*nrep*log(2*pi) - 0.5*log(detSigma) + log(likelihood)
  }
  logL
}

getMeanPostObsDeviance = function(mcmcfit, objects){
  mcmcLen = ((niter-nburnin)/nthin)
  
  nrep=objects[[1]]
  getRandMu = objects[[2]]
  ncomponents = objects[[3]]
  getEta = objects[[4]]
  nsubjects = objects[[5]]
  getSigma = objects[[6]]
  getObsLogL = objects[[7]]
  ds=objects[[8]]
  
  totalObsLogL=foreach(i=1:mcmcLen, .combine = c) %dopar%{
    
    sigma = getSigma(mcmcfit, mcmcIterNum = i)
    omega = solve(sigma)
    detSigma = det(sigma)
    
    randMu = getRandMu(mcmcfit, mcmcIterNum = i)
    betaTime = mcmcfit[[1]][i, "betaTime"]
    betaBy = mcmcfit[[1]][i, "betaBy"]
    betaGender = mcmcfit[[1]][i, "betaGender"]
    eta = getEta(mcmcfit, mcmcIterNum = i)
    
    -2*getObsLogL(nsubjects,nrep,ncomponents,ds$weight,
                  betaTime, ds$time, betaBy, as.numeric(ds$by)-1, betaGender,
                  as.numeric(ds$gender)-1, randMu, eta, 
                  detSigma, omega)
  }
  mean(totalObsLogL)
}

obsDIC=function(mcmcfit, summaryFuncName="mean"){
  meanPostDeviance = getMeanPostObsDeviance(mcmcfit,list(nrep, getRandMu, ncomponents,
                                                         getEta, nsubjects, getSigma,
                                                         getObsLogL, ds))

  eta = getEta(mcmcfit, summaryFuncName)
  randMu = getRandMu(mcmcfit, summaryFuncName)
  sigma = getSigma(mcmcfit, summaryFuncName)
  omega = solve(sigma)
  detSigma = det(sigma)
  
  randMu = getRandMu(mcmcfit, summaryFuncName)
  summaryFunc = get(summaryFuncName)
  betaTime = summaryFunc(mcmcfit[[1]][, "betaTime"])
  betaBy = summaryFunc(mcmcfit[[1]][, "betaBy"])
  betaGender = summaryFunc(mcmcfit[[1]][, "betaGender"])
  
  getPDandDIC(meanPostDeviance, -2*getObsLogL(nsubjects,nrep,ncomponents,ds$weight,
                                              betaTime, ds$time, betaBy, as.numeric(ds$by)-1, betaGender,
                                              as.numeric(ds$gender)-1, randMu, eta, 
                                              detSigma, omega))
}

calculateObsDIC1 = function(mcmcfit){
  obsDIC(mcmcfit, "mean")
}

calculateObsDIC2 = function(mcmcfit){
  obsDIC(mcmcfit, "modeFunc")
}

calculateObsDIC3 = function(mcmcfit, objects){
  
  nrep=objects[[1]]
  getRandMu = objects[[2]]
  ncomponents = objects[[3]]
  getEta = objects[[4]]
  nsubjects = objects[[5]]
  getSigma = objects[[6]]
  getObsLogL = objects[[7]]
  ds=objects[[8]]
  
  meanPostDeviance = getMeanPostObsDeviance(mcmcfit, objects)
  
  totalObsLogL=foreach(i=1:mcmcLen, .combine = c) %dopar%{
    
    sigma = getSigma(mcmcfit, mcmcIterNum = i)
    omega = solve(sigma)
    detSigma = det(sigma)
    
    randMu = getRandMu(mcmcfit, mcmcIterNum = i)
    betaTime = mcmcfit[[1]][i, "betaTime"]
    betaBy = mcmcfit[[1]][i, "betaBy"]
    betaGender = mcmcfit[[1]][i, "betaGender"]
    eta = getEta(mcmcfit, mcmcIterNum = i)
    
    getObsLogL(nsubjects,nrep,ncomponents,ds$weight,
                  betaTime, ds$time, betaBy, as.numeric(ds$by)-1, betaGender,
                  as.numeric(ds$gender)-1, randMu, eta, 
                  detSigma, omega)
  }
  #scale totalObsL for calculation purposes
  logfhaty = log(mean(exp(totalObsLogL)))
  print(exp(totalObsLogL))
  
  getPDandDIC(meanPostDeviance, -2*logfhaty)
}

calculateCompleteDataMeanPostDeviance=function(mcmcfit, objects){
  
  nrep=objects[[1]]
  getRandMu = objects[[2]]
  ncomponents = objects[[3]]
  getEta = objects[[4]]
  nsubjects = objects[[5]]
  getSigma = objects[[6]]
  getObsLogL = objects[[7]]
  ds=objects[[8]]
  
  dsweight = ds$weight
  dstime = ds$time
  dsgender = as.numeric(ds$gender)-1
  dsby = as.numeric(ds$by)-1
  
  mcmcLen = ((niter-nburnin)/nthin)
  
  totalDev=foreach(i=1:mcmcLen, .combine = c) %dopar%{
    ncomponents = objects[[3]] # doesn't work without keeping it here as well
    
    errPrecision = mcmcfit[[1]][i, "errPrecision"]
    randPrecision = mcmcfit[[1]][i, "randPrecision"]
    randMu = getRandMu(mcmcfit, mcmcIterNum = i)
    betaTime = mcmcfit[[1]][i, "betaTime"]
    betaBy = mcmcfit[[1]][i, "betaBy"]
    betaGender = mcmcfit[[1]][i, "betaGender"]
    eta = getEta(mcmcfit, mcmcIterNum = i)
    
    dev = 0
    for(j in 1:nsubjects){
      clusterNum = round(mcmcfit[[1]][i, paste("S[",j,"]", sep="")])
      randomIntercept = mcmcfit[[1]][i, paste("randomIntercept[",j,"]", sep="")]
      
      for(k in 1:nrep){
        index = (j-1)*nrep + k
        dev = dev + log(2*pi) - log(errPrecision)
        dev = dev + errPrecision * ((dsweight[index]-betaTime*dstime[index]-
                    betaBy*dsby[index]-betaGender*dsgender[index]-randomIntercept)^2)
        
        dev = dev - 2*log(eta[clusterNum])
        dev = dev + log(2*pi) - log(randPrecision)
        dev = dev + randPrecision * ((randomIntercept - randMu[clusterNum])^2)
      }
    }
    dev
  }
  mean(totalDev)
}

getAllocationMatrix=function(mcmcfit, mcmcLen, nsubjects){
  allocation = matrix(nrow=mcmcLen, ncol=nsubjects)
  for(i in 1:mcmcLen){
    for(j in 1:nsubjects){
      allocation[i,j] = round(mcmcfit[[1]][i, paste("S[",j,"]", sep="")])
    }
  }
  return(allocation)
}

getRandomInterceptMatrix=function(mcmcfit, mcmcLen, nsubjects){
  randomInterceptMatrix = matrix(nrow=mcmcLen, ncol=nsubjects)
  for(i in 1:mcmcLen){
    for(j in 1:nsubjects){
      randomInterceptMatrix[i, j] = mcmcfit[[1]][i, paste("randomIntercept[",j,"]", sep="")]
    }
  }
  return(randomInterceptMatrix)
}

calculateDIC4old=function(mcmcfit){
  meanPostDeviance = calculateCompleteDataMeanPostDeviance(mcmcfit,list(nrep, getRandMu, ncomponents,
                                                                        getEta, nsubjects, getSigma,
                                                                        getObsLogL, ds))
  mcmcLen = (niter-nburnin)/nthin
  allocation = getAllocationMatrix(mcmcfit, mcmcLen, nsubjects)
  randomInterceptMatrix = getRandomInterceptMatrix(mcmcfit, mcmcLen, nsubjects)
  
  secondHalfTheta=foreach(i=1:mcmcLen, .packages = 'R2jags') %dopar% {
    datanodes = list("dsweight"=ds$weight,"dsby"=as.numeric(ds$by)-1,
                     "dsgender"=as.numeric(ds$gender)-1,
                     "dstime"=ds$time,"nsubjects"=nsubjects,"nrep"=nrep, 
                  "varGammaParm"=0.0001,"betaTau"=0.0001, "betaMu"=0,
                  "ncomponents"=ncomponents, "dirichParm"=rep(1, ncomponents),
                  "S"=allocation[i,], "randomIntercept"= randomInterceptMatrix[i,])
    initialValues = list("betaGender"=c(0),"betaBy"=c(0), "betaTime"=c(0),
                         "errPrecision"=c(1),
                         "randPrecision"=rep(1, 1),
                         "Eta"=rep(1/ncomponents, ncomponents),
                         "mue"=quantile(extractRandomComp(viaReg = T), probs = seq(1/(ncomponents+1),ncomponents/(ncomponents+1), length.out = ncomponents)))
    stochasticNodes = c("betaGender","betaBy","betaTime", "errPrecision", 
                        "randPrecision", "randmu", "Eta")
    
    unload.module("glm")
    chainInits = list(initialValues)
    
    fit = jags(data=datanodes, inits=chainInits, stochasticNodes,
               n.chains=1, n.iter=niter, n.thin=nthin, n.burnin=nburnin,
               model.file=model, jags.module=NULL)
    mcmcfit = as.mcmc(fit)
  }
  
  dsweight = ds$weight
  dstime = ds$time
  dsgender = as.numeric(ds$gender)-1
  dsby = as.numeric(ds$by)-1
  
  totalDev=foreach(i=1:mcmcLen, .combine = c) %dopar%{
    mcmcLen = ((niter-nburnin)/nthin)
    
    mcmcfit = secondHalfTheta[[i]]
    
    errPrecision = modeFunc(mcmcfit[[1]][, "errPrecision"])
    randPrecision = modeFunc(mcmcfit[[1]][, "randPrecision"])
    randMu = getRandMu(mcmcfit, summaryFuncName = "modeFunc")
    betaTime = modeFunc(mcmcfit[[1]][, "betaTime"])
    betaBy = modeFunc(mcmcfit[[1]][, "betaBy"])
    betaGender = modeFunc(mcmcfit[[1]][, "betaGender"])
    eta = getEta(mcmcfit, summaryFuncName = "modeFunc")
    
    dev = 0
    for(j in 1:nsubjects){
      clusterNum = allocation[i,j]
      randomIntercept = randomInterceptMatrix[i,j]
      
      for(k in 1:nrep){
        index = (j-1)*nrep + k
        dev = dev + log(2*pi) - log(errPrecision)
        dev = dev + errPrecision * ((dsweight[index]-betaTime*dstime[index]-
                                       betaBy*dsby[index]-betaGender*dsgender[index]-randomIntercept)^2)
        
        dev = dev - 2*log(eta[clusterNum])
        dev = dev + log(2*pi) - log(randPrecision)
        dev = dev + randPrecision * ((randomIntercept - randMu[clusterNum])^2)
      }
    }
    dev
  }
  
  getPDandDIC(meanPostDeviance, mean(totalDev))
}


calculateDIC4=function(mcmcfit){
  meanPostDeviance = calculateCompleteDataMeanPostDeviance(mcmcfit,list(nrep, getRandMu, ncomponents,
                                                                        getEta, nsubjects, getSigma,
                                                                        getObsLogL, ds))
  mcmcLen = (niter-nburnin)/nthin
  allocation = getAllocationMatrix(mcmcfit, mcmcLen, nsubjects)
  randomInterceptMatrix = getRandomInterceptMatrix(mcmcfit, mcmcLen, nsubjects)
  
  secondHalfTheta=foreach(i=1:mcmcLen, .packages = 'R2jags',.combine = 'rbind') %dopar% {
    freqPerGroup = rep(0.000000001, ncomponents)
    tempFreq = table(allocation[i,])
    for(j in dimnames(tempFreq)[[1]]){
      freqPerGroup[as.numeric(j)] = tempFreq[j]
    }
    eta=(1+freqPerGroup)/(ncomponents + nsubjects)
    randMu = rep(0, ncomponents)
    for(j in 1:nsubjects){
      randMu[allocation[i,j]] = randMu[allocation[i,j]] + randomInterceptMatrix[i,j]
    }
    randMu = randMu/freqPerGroup
    
    withinGroupSsq = 0
    for(j in 1:nsubjects){
      withinGroupSsq = withinGroupSsq + (randomInterceptMatrix[i,j]-randMu[allocation[i,j]])^2
    }
    randPrecision = (nsubjects + 0.001 - 2)/(withinGroupSsq + 0.001)
    
    dsweight = ds$weight
    for(j in 1:nsubjects){
      for(k in 1:nrep){
        index = (j-1)*nrep + k
        dsweight[index] = dsweight[index] - randomInterceptMatrix[i,j]
      }
    }
    
    dstime = ds$time
    dsgender = as.numeric(ds$gender)-1
    dsby = as.numeric(ds$by)-1
    reg=lm(dsweight~dstime + dsgender + dsby + 0)
    errPrecision = (nsubjects*nrep-5)/sum((dsweight-reg$fitted.values)^2)
    
    list("eta"=eta, "randMu"=randMu, "randPrecision"=randPrecision, 
         "betaTime"=reg$coefficients["dstime"],"betaBy"=reg$coefficients["dsby"], 
         "betaGender"=reg$coefficients["dsgender"], "errPrecision"=errPrecision)
  }
  
  dsweight = ds$weight
  dstime = ds$time
  dsgender = as.numeric(ds$gender)-1
  dsby = as.numeric(ds$by)-1
  
  totalDev=foreach(i=1:mcmcLen, .combine = c) %dopar%{
    mcmcLen = ((niter-nburnin)/nthin)
    
    errPrecision = secondHalfTheta[i,]$errPrecision
    randPrecision = secondHalfTheta[i,]$randPrecision
    randMu = secondHalfTheta[i,]$randMu
    betaTime = secondHalfTheta[i,]$betaTime
    betaBy = secondHalfTheta[i,]$betaBy
    betaGender = secondHalfTheta[i,]$betaGender
    eta = secondHalfTheta[i,]$eta
    
    dev = 0
    for(j in 1:nsubjects){
      clusterNum = allocation[i,j]
      randomIntercept = randomInterceptMatrix[i,j]
      
      for(k in 1:nrep){
        index = (j-1)*nrep + k
        dev = dev + log(2*pi) - log(errPrecision)
        dev = dev + errPrecision * ((dsweight[index]-betaTime*dstime[index]-
                                       betaBy*dsby[index]-betaGender*dsgender[index]-randomIntercept)^2)
        
        dev = dev - 2*log(eta[clusterNum])
        dev = dev + log(2*pi) - log(randPrecision)
        dev = dev + randPrecision * ((randomIntercept - randMu[clusterNum])^2)
      }
    }
    dev
  }
  
  getPDandDIC(meanPostDeviance, mean(totalDev))
}


calculateDIC5 = function(mcmcfit){
  meanPostDeviance = calculateCompleteDataMeanPostDeviance(mcmcfit,list(nrep, getRandMu, ncomponents,
                                                                        getEta, nsubjects, getSigma,
                                                                        getObsLogL, ds))
  
  mcmcLen = (niter-nburnin)/nthin
  allocation = getAllocationMatrix(mcmcfit, mcmcLen, nsubjects)
  randomInterceptMatrix = getRandomInterceptMatrix(mcmcfit, mcmcLen, nsubjects)
  
  dsweight = ds$weight
  dstime = ds$time
  dsgender = as.numeric(ds$gender)-1
  dsby = as.numeric(ds$by)-1
  
  randMu=getRandMu(mcmcfit,summaryFuncName = "modeFunc")
  eta=getEta(mcmcfit,summaryFuncName = "modeFunc") 
  errPrecision = modeFunc(mcmcfit[[1]][, "errPrecision"])
  randPrecision = modeFunc(mcmcfit[[1]][, "randPrecision"])
  betaTime = modeFunc(mcmcfit[[1]][, "betaTime"])
  betaBy = modeFunc(mcmcfit[[1]][, "betaBy"])
  betaGender = modeFunc(mcmcfit[[1]][, "betaGender"])
  
  totalDev=foreach(j=1:nsubjects, .combine = c) %dopar%{
    dev = 0
    clusterNum = round(modeFunc(allocation[,j]))
    randomIntercept = modeFunc(randomInterceptMatrix[,j])
    
    for(k in 1:nrep){
      index = (j-1)*nrep + k
      dev = dev + log(2*pi) - log(errPrecision)
      dev = dev + errPrecision * ((dsweight[index]-betaTime*dstime[index]-
                                     betaBy*dsby[index]-betaGender*dsgender[index]-randomIntercept)^2)
      
      dev = dev - 2*log(eta[clusterNum])
      dev = dev + log(2*pi) - log(randPrecision)
      dev = dev + randPrecision * ((randomIntercept - randMu[clusterNum])^2)
    }
    dev
  }
  
  getPDandDIC(meanPostDeviance, sum(totalDev))
}

#following function is based on f(y|b,s), MAP means maximum a posteriori
calculateDIC7 = function(mcmcfit){
  
  errMAP=modeFunc(mcmcfit[[1]][, paste("errPrecision", sep="")])
  betaTime = modeFunc(mcmcfit[[1]][, "betaTime"])
  betaBy = modeFunc(mcmcfit[[1]][, "betaBy"])
  betaGender = modeFunc(mcmcfit[[1]][, "betaGender"])
  
  dsweight = ds$weight
  dstime = ds$time
  dsgender = as.numeric(ds$gender)-1
  dsby = as.numeric(ds$by)-1
  
  totalDev=foreach(i=1:nsubjects, .combine = c) %dopar%{
    dev = 0
    randomIntMAP=modeFunc(mcmcfit[[1]][, paste("randomIntercept[",i,"]", sep="")])
    for(j in 1:nrep){
      index = (i-1)*nrep + j
      bias = dsweight[index] - betaTime*dstime[index] - betaBy*dsby[index] - 
        betaGender*dsgender[index] - randomIntMAP
      dev = dev + log(2*pi) - log(errMAP) + errMAP*(bias^2)
    }
    dev
  }
  
  Dthetabar = sum(totalDev)
  
  mcmcLen = ((niter-nburnin)/nthin)
  totalDev=foreach(k=1:mcmcLen, .combine = c) %dopar%{
    errMAP=mcmcfit[[1]][k, paste("errPrecision", sep="")]
    betaTime = mcmcfit[[1]][k, "betaTime"]
    betaBy = mcmcfit[[1]][k, "betaBy"]
    betaGender = mcmcfit[[1]][k, "betaGender"]
    
    dev=0
    for(i in 1:nsubjects){
      randomIntMAP=mcmcfit[[1]][k, paste("randomIntercept[",i,"]", sep="")]
      for(j in 1:nrep){
        index = (i-1)*nrep + j
        bias = dsweight[index] - betaTime*dstime[index] - betaBy*dsby[index] - 
          betaGender*dsgender[index] - randomIntMAP
        dev = dev + log(2*pi) - log(errMAP) + errMAP*(bias^2)
      }
    }
    dev
  }
  
  getPDandDIC(mean(totalDev), Dthetabar)
}
