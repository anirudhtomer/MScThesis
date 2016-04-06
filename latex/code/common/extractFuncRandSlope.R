modeFunc=function(data){
  d = density(data, n=1e4)
  return(d$x[which.max(d$y)])
}

getEta=function(mcmcfit, summaryFuncName="mean", mcmcIterNum = NA){
  ncomponents = attributes(mcmcfit)$ncomponents
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
  ncomponents = attributes(mcmcfit)$ncomponents
  summaryFunc = get(summaryFuncName)
  randmu=list()
  
  for(i in 1:ncomponents){
    if(!is.numeric(mcmcIterNum)){
      randmu[[i]] = c(summaryFunc(mcmcfit[[1]][,paste("randmu[",i,",1]", sep="")]),
                      summaryFunc(mcmcfit[[1]][,paste("randmu[",i,",2]", sep="")]))
    }else{
      randmu[[i]] = c(mcmcfit[[1]][mcmcIterNum, paste("randmu[",i,",1]", sep="")],
                      mcmcfit[[1]][mcmcIterNum, paste("randmu[",i,",2]", sep="")])
    }
  }
  return(randmu)
}

getSigma=function(mcmcfit,  summaryFuncName="mean", mcmcIterNum = NA){
  ncomponents = attributes(mcmcfit)$ncomponents
  nrep = attributes(mcmcfit)$nrep
  
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
  ncomponents = attributes(mcmcfit)$ncomponents
  nrep = attributes(mcmcfit)$nrep
  
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

getAllocations=function(mcmcfit, summaryFuncName = "mean", mcmcIterNum=NA){
  nsubjects = attributes(mcmcfit)$nsubjects
  allocation = numeric(nsubjects)
  summaryFunc = get(summaryFuncName)
  
  if(!is.numeric(mcmcIterNum)){
    for(j in 1:nsubjects){
      allocation[j] = round(summaryFunc(mcmcfit[[1]][, paste("S[",j,"]", sep="")]))
    }
  }else{
    for(j in 1:nsubjects){
      allocation[j] = mcmcfit[[1]][mcmcIterNum, paste("S[",j,"]", sep="")]
    }
  }
  return(allocation)
}

getXBeta=function(mcmcfit,  summaryFuncName="mean", mcmcIterNum = NA){
  nsubjects = attributes(mcmcfit)$nsubjects
  nrep = attributes(mcmcfit)$nrep
  
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