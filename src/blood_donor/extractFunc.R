source("../common/extractFuncRandSlopeConditional.R")

getSigma=function(mcmcfit,  summaryFuncName="mean", donLast2Yrs, mcmcIterNum = NA){
  ncomponents = attributes(mcmcfit)$ncomponents
  nrep = length(donLast2Yrs)
  time=attributes(mcmcfit)$time
  INTERCEPT_SCALE=attributes(mcmcfit)$INTERCEPT_SCALE
  
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
        sigma[[b]][i,j] = d22*donLast2Yrs[i]*donLast2Yrs[j] + INTERCEPT_SCALE*(donLast2Yrs[i]+donLast2Yrs[j])*d12 + d11*INTERCEPT_SCALE*INTERCEPT_SCALE
        if(i==j){
          sigma[[b]][i,j] = sigma[[b]][i,j] + errVariance
        }
      }
    }
  }
  
  return(sigma)
}

extractRandomComp = function(){
  temp = ds$Hb
  reg=lm(Hb~Age + Season + Donate + TSPD +  donationLast2Years + 
           donateLast2Donate + donateLast2TSPD + donateLast2Square, data=ds)
  
  temp = reg$residuals
  
  randIntercept = numeric()
  randSlope = numeric()
  prevEnd = 0
  k=0
  for(i in unique(ds$Id)){
    startIndex = prevEnd + 1
    endIndex = prevEnd + length(timePerSubject[[paste(i)]])
    
    reg=lm(temp[startIndex:endIndex]~0+I(rep(0.1,endIndex-startIndex+1)) + I(donationLast2YearsPerSubject[[paste(i)]]))   
    randIntercept[k] = reg$coefficients[1]
    randSlope[k] = reg$coefficients[2]
    k=k+1
    
    prevEnd = endIndex
  }
  return(cbind(randIntercept, randSlope))
}

randomComp=extractRandomComp()
randomCompDf = data.frame(randomComp[complete.cases(randomComp)==TRUE,])
qplot(x=randIntercept, y=randSlope, data=randomCompDf, xlab="Random intercept", ylab="Random slope")
#var(randomCompDf)


