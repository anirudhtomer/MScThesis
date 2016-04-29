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
  sigma = list(ncomponents)
  
  for(b in 1:ncomponents){
    sigma[[b]]=matrix(nrow=nrep, ncol=nrep)
    if(!is.numeric(mcmcIterNum)){
      for(i in 1:nrep){
        for(j in 1:nrep){
          sigma[[b]][i,j] = summaryFunc(mcmcfit[[1]][, paste("sigma[",i,",",j,",",b,"]", sep="")])
        }
      }
    }else{
      for(i in 1:nrep){
        for(j in 1:nrep){
          sigma[[b]][i,j] = mcmcfit[[1]][mcmcIterNum, paste("sigma[",i,",",j,",",b,"]", sep="")]
        }
      }
    }
  }
  return(sigma)
}

getOmega=function(mcmcfit,  summaryFuncName="mean", mcmcIterNum = NA){
  summaryFunc = get(summaryFuncName)
  omega = list(ncomponents)
  
  for(b in 1:ncomponents){
    omega[[b]]=matrix(nrow=nrep, ncol=nrep)
    if(!is.numeric(mcmcIterNum)){
      for(i in 1:nrep){
        for(j in 1:nrep){
          omega[[b]][i,j] = summaryFunc(mcmcfit[[1]][, paste("omega[",i,",",j,",",b,"]", sep="")])
        }
      }
    }else{
      for(i in 1:nrep){
        for(j in 1:nrep){
          omega[[b]][i,j] = mcmcfit[[1]][mcmcIterNum, paste("omega[",i,",",j,",",b,"]", sep="")]
        }
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

mcmcLen = nrow(mcmcfit[[1]])
detCov=foreach(i=1:mcmcLen, .combine='rbind', .packages='MASS') %dopar%{
  eta = getEta(mcmcfit, mcmcIterNum = i)
  randmu = getRandMu(mcmcfit, mcmcIterNum = i)
  sigmaEst = getSigma(mcmcfit, mcmcIterNum = i)
  xBeta = getXBeta(mcmcfit, mcmcIterNum = i)
  
  sampleRandomPart = dswide_y-xBeta
  
  newRandomPart = matrix(nrow=0, ncol=nrep)
  for(j in 1:ncomponents){
    sampleSize = floor(eta[j]*180)
    if(sampleSize>0){
      newRandomPart=rbind(newRandomPart, mvrnorm(n=sampleSize, mu = rep(randmu[j], nrep), Sigma = sigmaEst[[j]]))
    }
  }
  
  newRange = range(newRandomPart)
  sampleRange = range(sampleRandomPart)
  
  #return(matrix(c(mean(newRandomPart), mean(sampleRandomPart), T,F), nrow=2))
  return(matrix(c(newRange[2]-newRange[1], sampleRange[2]-sampleRange[1], T,F), nrow=2))
}
detCovDf = data.frame(detCov)
detCovDf$X2 = factor(detCovDf$X2,labels = c("sample", "new"))
HPDinterval(mcmc(detCovDf[detCovDf$X2=="sample",]$X1))
qplot(X1, data = detCovDf, geom = "density", group=X2, color=X2)
beep(sound = 8)

sum(detCov[(1:2000)%%2==0,1] > detCov[(1:2000)%%2==1,1])


mahalaDist=foreach(i=1:mcmcLen, .combine='rbind') %do%{
  eta = getEta(mcmcfit, mcmcIterNum = i)
  randmu = getRandMu(mcmcfit, mcmcIterNum = i)
  sigmaEst = getSigma(mcmcfit, mcmcIterNum = i)
  xBeta = getXBeta(mcmcfit, mcmcIterNum = i)
  
  obsMahal = rep(0, nsubjects)
  #for(j in 1:nsubjects){
   # alloc = mcmcfit[[1]][i, paste("S[",j, "]", sep="")]
    #randomPart = dswide_y[j,]-xBeta[j,]
    #for(m in 1:ncomponents){
    #  if(alloc!=m){
     #   obsMahal[j] = obsMahal[j] + mahalanobis(randomPart, rep(randmu[m], nrep), 
                                               # sigmaEst)
    #  }
    #}
    #obsMahal[j] = obsMahal[j]/(ncomponents-1)
  #}
  
  newMahal = numeric(0)
  newObsCount = 0
  for(alloc in 1:ncomponents){
    newObs=mvrnorm(n=round(eta[alloc]*1000), rep(randmu[alloc], nrep), sigmaEst)
    for(f in 1:nrow(newObs)){
      newObsCount = newObsCount + 1
      newMahal[newObsCount] = 0
      for(g in 1:ncomponents){
        if(g!=alloc){
          newMahal[newObsCount] = newMahal[newObsCount] + 
            mahalanobis(newObs[f,], rep(randmu[g], nrep), sigmaEst)
        }
      }
      newMahal[newObsCount] = newMahal[newObsCount]/(ncomponents-1)
    }
  }

  return(mean(newMahal))
}
beep(sound=8)
HPDinterval(mcmc(mahalaDist))
quantile(mahalaDist,probs = c(0.025, 0.975))
plot(density(mahalaDist))
