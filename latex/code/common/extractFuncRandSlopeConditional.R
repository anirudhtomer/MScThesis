source("../common/extractFuncRandSlope.R")

getRandomComp=function(mcmcfit, summaryFuncName="mean", mcmcIterNum=NA){
  summaryFunc = get(summaryFuncName)
  
  nsubjects = attributes(mcmcfit)$nsubjects
  randomComp = matrix(nrow = nsubjects, ncol=2)
  
  if(!is.numeric(mcmcIterNum)){
    for(j in 1:nsubjects){
      randomComp[j,1] = summaryFunc(mcmcfit[[1]][, paste("randomComp[",j,",1]", sep="")])
      randomComp[j,2] = summaryFunc(mcmcfit[[1]][, paste("randomComp[",j,",2]", sep="")])
    }
  }else{
    for(j in 1:nsubjects){
      randomComp[j,1] = mcmcfit[[1]][mcmcIterNum, paste("randomComp[",j,",1]", sep="")]
      randomComp[j,2] = mcmcfit[[1]][mcmcIterNum, paste("randomComp[",j,",2]", sep="")]
    }
  }
  return(randomComp)
}

getRandSigma=function(mcmcfit, summaryFuncName="mean", mcmcIterNum=NA){
  ncomponents = attributes(mcmcfit)$ncomponents
  
  summaryFunc = get(summaryFuncName)
  randSigma = list(ncomponents)
  
  for(b in 1:ncomponents){
    d11=d22=d12=0
    
    if(!is.numeric(mcmcIterNum)){
      d11 = summaryFunc(mcmcfit[[1]][,paste("randSigma[1,1,",b,"]", sep="")])
      d22 = summaryFunc(mcmcfit[[1]][,paste("randSigma[2,2,",b,"]", sep="")])
      d12 = summaryFunc(mcmcfit[[1]][,paste("randSigma[1,2,",b,"]", sep="")])
    }else{
      d11 = mcmcfit[[1]][mcmcIterNum, paste("randSigma[1,1,",b,"]", sep="")]
      d22 = mcmcfit[[1]][mcmcIterNum, paste("randSigma[2,2,",b,"]", sep="")]
      d12 = mcmcfit[[1]][mcmcIterNum, paste("randSigma[1,2,",b,"]", sep="")]
    }
    randSigma[[b]]=matrix(c(d11,d12,d12,d22),nrow=2, ncol=2)
  }
  
  return(randSigma)
}

getRandPrecision=function(mcmcfit, summaryFuncName="mean", mcmcIterNum=NA){
  ncomponents = attributes(mcmcfit)$ncomponents
  
  summaryFunc = get(summaryFuncName)
  randPrecision = list(ncomponents)
  
  for(b in 1:ncomponents){
    d11=d22=d12=0
    
    if(!is.numeric(mcmcIterNum)){
      d11 = summaryFunc(mcmcfit[[1]][,paste("randPrecision[1,1,",b,"]", sep="")])
      d22 = summaryFunc(mcmcfit[[1]][,paste("randPrecision[2,2,",b,"]", sep="")])
      d12 = summaryFunc(mcmcfit[[1]][,paste("randPrecision[1,2,",b,"]", sep="")])
    }else{
      d11 = mcmcfit[[1]][mcmcIterNum, paste("randPrecision[1,1,",b,"]", sep="")]
      d22 = mcmcfit[[1]][mcmcIterNum, paste("randPrecision[2,2,",b,"]", sep="")]
      d12 = mcmcfit[[1]][mcmcIterNum, paste("randPrecision[1,2,",b,"]", sep="")]
    }
    randPrecision[[b]]=matrix(c(d11,d12,d12,d22),nrow=2, ncol=2)
  }
  
  return(randPrecision)
}

getSigma=function(mcmcfit,  summaryFuncName="mean", mcmcIterNum = NA){
  ncomponents = attributes(mcmcfit)$ncomponents
  nrep = attributes(mcmcfit)$nrep
  time=attributes(mcmcfit)$time
  
  summaryFunc = get(summaryFuncName)
  sigma = list(ncomponents)
  
  for(b in 1:ncomponents){
    sigma[[b]]=matrix(nrow=nrep, ncol=nrep)

    errVariance=d11=d22=d12=0
    
    if(!is.numeric(mcmcIterNum)){
      errVariance = summaryFunc(1/mcmcfit[[1]][,"errPrecision"])
      d11 = summaryFunc(mcmcfit[[1]][,paste("randSigma[1,1,",b,"]", sep="")])
      d22 = summaryFunc(mcmcfit[[1]][,paste("randSigma[2,2,",b,"]", sep="")])
      d12 = summaryFunc(mcmcfit[[1]][,paste("randSigma[1,2,",b,"]", sep="")])
    }else{
      errVariance = 1/mcmcfit[[1]][mcmcIterNum,"errPrecision"]
      d11 = mcmcfit[[1]][mcmcIterNum, paste("randSigma[1,1,",b,"]", sep="")]
      d22 = mcmcfit[[1]][mcmcIterNum, paste("randSigma[2,2,",b,"]", sep="")]
      d12 = mcmcfit[[1]][mcmcIterNum, paste("randSigma[1,2,",b,"]", sep="")]
    }

    for(i in 1:nrep){
      for(j in 1:nrep){
        sigma[[b]][i,j] = d22*time[i]*time[j] + (time[i]+time[j])*d12 + d11
        if(i==j){
          sigma[[b]][i,j] = sigma[[b]][i,j] + errVariance
        }
      }
    }
  }
  
  return(sigma)
}