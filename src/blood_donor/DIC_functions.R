library(mvtnorm)
 
INTERCEPT_SCALE = 0.1

getPDandDIC=function(meanPostDeviance, Dthetabar){
  pd = meanPostDeviance - Dthetabar
  dic = meanPostDeviance + pd
  
  list("pd"=pd,"dic"=dic)
}

calculateObsDThetaBar=function(mcmcfit, summaryFuncName){
  summaryFunc = get(summaryFuncName)
  
  ncomponents = attributes(mcmcfit)$ncomponents
  nsubjects = attributes(mcmcfit)$nsubjects
  
  etaSummary = getEta(mcmcfit, summaryFuncName)
  randmuSummary = getRandMu(mcmcfit, summaryFuncName)
  
  betaAge = summaryFunc(mcmcfit[[1]][,"betaAge"])
  betaDonate = summaryFunc(mcmcfit[[1]][,"betaDonate"])
  betaSeason = summaryFunc(mcmcfit[[1]][,"betaSeason"])
  betaTSPD = summaryFunc(mcmcfit[[1]][,"betaTSPD"])
  betaDonateLast2TSPD = summaryFunc(mcmcfit[[1]][,"betaDonateLast2TSPD"])
  betaDonateLast2Donate = summaryFunc(mcmcfit[[1]][,"betaDonateLast2Donate"])
  betaDonateLast2Square = summaryFunc(mcmcfit[[1]][,"betaDonateLast2Square"])
  
  Dthetabar = 0
  for(i in 1:nsubjects){
    temp = 0
    startIndex = cumsumHb[i]+1
    endIndex = cumsumHb[i]+numHb[i]
    
    sigmaSummary = getSigma(mcmcfit, summaryFuncName, donLast2Yrs = ds$donationLast2Years[startIndex:endIndex])
    fixedPartMean = betaAge*ds$Age[startIndex:endIndex] + betaSeason*(as.numeric(ds$Season[startIndex:endIndex])-1) + 
      betaDonate*ds$Donate[startIndex:endIndex] + betaTSPD*ds$TSPD[startIndex:endIndex] + 
      betaDonateLast2TSPD*ds$donateLast2TSPD[startIndex:endIndex] +
      betaDonateLast2Donate*ds$donateLast2Donate[startIndex:endIndex] +
      betaDonateLast2Square*ds$donateLast2Square[startIndex:endIndex]
    
    for(j in 1:ncomponents){
      randomPartMean = randmuSummary[[j]][1]*INTERCEPT_SCALE + randmuSummary[[j]][2]*ds$donationLast2Years[startIndex:endIndex]
      temp = temp + etaSummary[j] * dmvnorm(x = ds$Hb[startIndex:endIndex], 
                                            mean = randomPartMean + fixedPartMean, 
                                            sigma = sigmaSummary[[j]], log = F)
    }
    Dthetabar = Dthetabar + log(temp)
  }
  Dthetabar = -2*Dthetabar
}

calculateObsTotalDTheta = function(mcmcfit){
  ds = ds
  cumsumHb = cumsumHb
  numHb = numHb
  INTERCEPT_SCALE = INTERCEPT_SCALE
  
  ncomponents = attributes(mcmcfit)$ncomponents
  nsubjects = attributes(mcmcfit)$nsubjects
  
  totalDtheta=foreach(k=1:nrow(mcmcfit[[1]]), .combine = c, .packages = c('mvtnorm')) %dopar%{
    
    source("extractFunc.R")
    
    eta = getEta(mcmcfit, mcmcIterNum = k)
    randmu = getRandMu(mcmcfit, mcmcIterNum = k)
    
    betaAge = mcmcfit[[1]][k,"betaAge"]
    betaDonate = mcmcfit[[1]][k,"betaDonate"]
    betaSeason = mcmcfit[[1]][k,"betaSeason"]
    betaTSPD = mcmcfit[[1]][k,"betaTSPD"]
    betaDonateLast2TSPD = mcmcfit[[1]][k,"betaDonateLast2TSPD"]
    betaDonateLast2Donate = mcmcfit[[1]][k,"betaDonateLast2Donate"]
    betaDonateLast2Square = mcmcfit[[1]][k,"betaDonateLast2Square"]
    
    Dtheta = 0
    for(i in 1:nsubjects){
      temp = 0
      startIndex = cumsumHb[i]+1
      endIndex = cumsumHb[i]+numHb[i]
      fixedPartMean = betaAge*ds$Age[startIndex:endIndex] + betaSeason*(as.numeric(ds$Season[startIndex:endIndex])-1) + 
        betaDonate*ds$Donate[startIndex:endIndex] + betaTSPD*ds$TSPD[startIndex:endIndex] + 
        betaDonateLast2TSPD*ds$donateLast2TSPD[startIndex:endIndex] +
        betaDonateLast2Donate*ds$donateLast2Donate[startIndex:endIndex] +
        betaDonateLast2Square*ds$donateLast2Square[startIndex:endIndex]
      
      sigma = getSigma(mcmcfit, donLast2Yrs = ds$donationLast2Years[startIndex:endIndex], mcmcIterNum = k)
      for(j in 1:ncomponents){
        randomPartMean = randmu[[j]][1]*INTERCEPT_SCALE + randmu[[j]][2]*ds$donationLast2Years[startIndex:endIndex]
        
        temp = temp + eta[j] * dmvnorm(x = ds$Hb[startIndex:endIndex], 
                                       mean = randomPartMean + fixedPartMean, 
                                       sigma = sigma[[j]], log = F)
      }
      Dtheta = Dtheta + log(temp)
    }
    -2*Dtheta
  }
}

calculateDIC1 = function(mcmcfit){
  getPDandDIC(mean(calculateObsTotalDTheta(mcmcfit)), calculateObsDThetaBar(mcmcfit, "mean"))
}

calculateDIC2 = function(mcmcfit){
  getPDandDIC(mean(calculateObsTotalDTheta(mcmcfit)), calculateObsDThetaBar(mcmcfit, "modeFunc"))
}

calculateDIC3 = function(mcmcfit){
  
  DThetabar_functionalapprox = 0
  
  ds = ds
  cumsumHb = cumsumHb
  numHb = numHb
  INTERCEPT_SCALE = INTERCEPT_SCALE
  
  ncomponents = attributes(mcmcfit)$ncomponents
  nsubjects = attributes(mcmcfit)$nsubjects

  obsDev_Subjects=foreach(k=1:nrow(mcmcfit[[1]]), .combine = 'rbind', .packages = c('mvtnorm')) %dopar%{
    
    source("extractFunc.R")
    
    eta = getEta(mcmcfit, mcmcIterNum = k)
    randmu = getRandMu(mcmcfit, mcmcIterNum = k)
    
    betaAge = mcmcfit[[1]][k,"betaAge"]
    betaDonate = mcmcfit[[1]][k,"betaDonate"]
    betaSeason = mcmcfit[[1]][k,"betaSeason"]
    betaTSPD = mcmcfit[[1]][k,"betaTSPD"]
    betaDonateLast2TSPD = mcmcfit[[1]][k,"betaDonateLast2TSPD"]
    betaDonateLast2Donate = mcmcfit[[1]][k,"betaDonateLast2Donate"]
    betaDonateLast2Square = mcmcfit[[1]][k,"betaDonateLast2Square"]
    
    Dtheta = rep(0, nsubjects)
    for(i in 1:nsubjects){
      startIndex = cumsumHb[i]+1
      endIndex = cumsumHb[i]+numHb[i]
      
      fixedPartMean = betaAge*ds$Age[startIndex:endIndex] + betaSeason*(as.numeric(ds$Season[startIndex:endIndex])-1) + 
        betaDonate*ds$Donate[startIndex:endIndex] + betaTSPD*ds$TSPD[startIndex:endIndex] + 
        betaDonateLast2TSPD*ds$donateLast2TSPD[startIndex:endIndex] +
        betaDonateLast2Donate*ds$donateLast2Donate[startIndex:endIndex] +
        betaDonateLast2Square*ds$donateLast2Square[startIndex:endIndex]
      
      sigma = getSigma(mcmcfit, donLast2Yrs = ds$donationLast2Years[startIndex:endIndex], mcmcIterNum = k)
      for(j in 1:ncomponents){
        randomPartMean = randmu[[j]][1]*INTERCEPT_SCALE + randmu[[j]][2]*ds$donationLast2Years[startIndex:endIndex]
       
        Dtheta[i] = Dtheta[i] + eta[j] * dmvnorm(x = ds$Hb[startIndex:endIndex], 
                                           mean = randomPartMean + fixedPartMean, 
                                           sigma = sigma[[j]], log = F)
      }
    }
    Dtheta
  }
  
  DThetabar_functionalapprox = -2*sum(log(apply(X = obsDev_Subjects, MARGIN = 2, FUN = mean)))
  getPDandDIC(mean(calculateObsTotalDTheta(mcmcfit)), DThetabar_functionalapprox)
}

calculateDIC7 = function(mcmcfit, summaryFuncName="modeFunc"){
  meanPostDeviance = 0
  
  ds = ds
  cumsumHb = cumsumHb
  numHb = numHb
  INTERCEPT_SCALE = INTERCEPT_SCALE
  
  ncomponents = attributes(mcmcfit)$ncomponents
  nsubjects = attributes(mcmcfit)$nsubjects

  totalDtheta=foreach(k=1:nrow(mcmcfit[[1]]), .combine = 'c', .packages = c('mvtnorm')) %dopar%{
    source("extractFunc.R")
    
    betaAge = mcmcfit[[1]][k,"betaAge"]
    betaDonate = mcmcfit[[1]][k,"betaDonate"]
    betaSeason = mcmcfit[[1]][k,"betaSeason"]
    betaTSPD = mcmcfit[[1]][k,"betaTSPD"]
    betaDonateLast2TSPD = mcmcfit[[1]][k,"betaDonateLast2TSPD"]
    betaDonateLast2Donate = mcmcfit[[1]][k,"betaDonateLast2Donate"]
    betaDonateLast2Square = mcmcfit[[1]][k,"betaDonateLast2Square"]
    
    randComp = getRandomComp(mcmcfit, mcmcIterNum = k)
    errVariance = 1/mcmcfit[[1]][k,"errPrecision"]
    
    Dtheta = 0
    for(i in 1:nsubjects){
      startIndex = cumsumHb[i]+1
      endIndex = cumsumHb[i]+numHb[i]
      
      fixedPartMean = betaAge*ds$Age[startIndex:endIndex] + betaSeason*(as.numeric(ds$Season[startIndex:endIndex])-1) + 
        betaDonate*ds$Donate[startIndex:endIndex] + betaTSPD*ds$TSPD[startIndex:endIndex] + 
        betaDonateLast2TSPD*ds$donateLast2TSPD[startIndex:endIndex] +
        betaDonateLast2Donate*ds$donateLast2Donate[startIndex:endIndex] +
        betaDonateLast2Square*ds$donateLast2Square[startIndex:endIndex]
      
      randomPartMean = randComp[i,1]*INTERCEPT_SCALE + randComp[i,2]*ds$donationLast2Years[startIndex:endIndex]
      Dtheta = Dtheta + dmvnorm(x = ds$Hb[startIndex:endIndex], 
                                     mean = randomPartMean + fixedPartMean, 
                                     sigma = diag(numHb[i])*errVariance, log = T)
    }
    Dtheta
  }
  meanPostDeviance = -2*mean(totalDtheta)
  
  #Calcuing D(thetabar)
  summaryFunc = get(summaryFuncName)
  
  randComp = getRandomComp(mcmcfit, summaryFuncName)
  errVariance = summaryFunc(1/mcmcfit[[1]][,"errPrecision"])
  betaAge = summaryFunc(mcmcfit[[1]][,"betaAge"])
  betaDonate = summaryFunc(mcmcfit[[1]][,"betaDonate"])
  betaSeason = summaryFunc(mcmcfit[[1]][,"betaSeason"])
  betaTSPD = summaryFunc(mcmcfit[[1]][,"betaTSPD"])
  betaDonateLast2TSPD = summaryFunc(mcmcfit[[1]][,"betaDonateLast2TSPD"])
  betaDonateLast2Donate = summaryFunc(mcmcfit[[1]][,"betaDonateLast2Donate"])
  betaDonateLast2Square = summaryFunc(mcmcfit[[1]][,"betaDonateLast2Square"])
  
  Dthetabar = 0
  for(i in 1:nsubjects){
    startIndex = cumsumHb[i]+1
    endIndex = cumsumHb[i]+numHb[i]
    
    fixedPartMean = betaAge*ds$Age[startIndex:endIndex] + betaSeason*(as.numeric(ds$Season[startIndex:endIndex])-1) + 
      betaDonate*ds$Donate[startIndex:endIndex] + betaTSPD*ds$TSPD[startIndex:endIndex] + 
      betaDonateLast2TSPD*ds$donateLast2TSPD[startIndex:endIndex] +
      betaDonateLast2Donate*ds$donateLast2Donate[startIndex:endIndex] +
      betaDonateLast2Square*ds$donateLast2Square[startIndex:endIndex]
    
    randomPartMean = randComp[i,1]*INTERCEPT_SCALE + randComp[i,2]*ds$donationLast2Years[startIndex:endIndex]
    Dthetabar = Dthetabar + dmvnorm(x = ds$Hb[startIndex:endIndex], 
                              mean = randomPartMean + fixedPartMean, 
                              sigma = diag(numHb[i])*errVariance, log = T)
  }
  
  getPDandDIC(meanPostDeviance, -2*Dthetabar)
}

calculateDIC5 = function(mcmcfit){
  summaryFuncName = "modeFunc"
  
  #Calcuing D(thetabar)
  summaryFunc = get(summaryFuncName)
  ncomponents = attributes(mcmcfit)$ncomponents
  nsubjects = attributes(mcmcfit)$nsubjects
  INTERCEPT_SCALE = INTERCEPT_SCALE

  eta = getEta(mcmcfit, summaryFuncName)
  randComp = getRandomComp(mcmcfit, summaryFuncName)
  randmu = getRandMu(mcmcfit, summaryFuncName)
  randSigma = getRandSigma(mcmcfit, summaryFuncName)
  
  betaAge = summaryFunc(mcmcfit[[1]][,"betaAge"])
  betaDonate = summaryFunc(mcmcfit[[1]][,"betaDonate"])
  betaSeason = summaryFunc(mcmcfit[[1]][,"betaSeason"])
  betaTSPD = summaryFunc(mcmcfit[[1]][,"betaTSPD"])
  betaDonateLast2TSPD = summaryFunc(mcmcfit[[1]][,"betaDonateLast2TSPD"])
  betaDonateLast2Donate = summaryFunc(mcmcfit[[1]][,"betaDonateLast2Donate"])
  betaDonateLast2Square = summaryFunc(mcmcfit[[1]][,"betaDonateLast2Square"])
  
  errVariance = summaryFunc(1/mcmcfit[[1]][,"errPrecision"])
  allocations = getAllocations(mcmcfit, summaryFuncName)
  
  Dthetabar = 0
  for(i in 1:nsubjects){
    startIndex = cumsumHb[i]+1
    endIndex = cumsumHb[i]+numHb[i]
    
    randomPartMean = randComp[i,1]*INTERCEPT_SCALE + randComp[i,2]*ds$donationLast2Years[startIndex:endIndex]
    
    fixedPartMean = betaAge*ds$Age[startIndex:endIndex] + betaSeason*(as.numeric(ds$Season[startIndex:endIndex])-1) + 
      betaDonate*ds$Donate[startIndex:endIndex] + betaTSPD*ds$TSPD[startIndex:endIndex] + 
      betaDonateLast2TSPD*ds$donateLast2TSPD[startIndex:endIndex] +
      betaDonateLast2Donate*ds$donateLast2Donate[startIndex:endIndex] +
      betaDonateLast2Square*ds$donateLast2Square[startIndex:endIndex]
    
    Dthetabar = Dthetabar + log(eta[allocations[i]]) +
      dmvnorm(x = ds$Hb[startIndex:endIndex],
              mean = randomPartMean + fixedPartMean, 
              sigma = diag(numHb[i])*errVariance, log = T) +
      dmvnorm(x = randComp[i,],
              mean = randmu[[allocations[i]]], 
              sigma = randSigma[[allocations[i]]], log = T)
  }
  
  getPDandDIC(calculateCompleteMeanPostDeviance(mcmcfit), -2*Dthetabar)
}

calculateDIC4 = function(mcmcfit){
  E_theta_givenY_Z=calculate_E_theta_givenY_Z(mcmcfit)
  
  ncomponents = attributes(mcmcfit)$ncomponents
  nsubjects = attributes(mcmcfit)$nsubjects

  ds = ds
  cumsumHb = cumsumHb
  numHb = numHb
  INTERCEPT_SCALE=INTERCEPT_SCALE
  
  #Calcuing D(thetabar)
  totalDthetabar = foreach(k=1:nrow(mcmcfit[[1]]), .combine = 'c', .packages = c('mvtnorm')) %dopar%{
    eta = E_theta_givenY_Z[[k]]$eta
    randmu = E_theta_givenY_Z[[k]]$randmu
    randSigma = E_theta_givenY_Z[[k]]$randSigma
    
    betaAge = E_theta_givenY_Z[[k]]$betaAge
    betaDonate = E_theta_givenY_Z[[k]]$betaDonate
    betaSeason = E_theta_givenY_Z[[k]]$betaSeason
    betaTSPD = E_theta_givenY_Z[[k]]$betaTSPD
    betaDonateLast2TSPD = E_theta_givenY_Z[[k]]$betaDonateLast2TSPD
    betaDonateLast2Donate = E_theta_givenY_Z[[k]]$betaDonateLast2Donate
    betaDonateLast2Square = E_theta_givenY_Z[[k]]$betaDonateLast2Square
 
    errVariance = E_theta_givenY_Z[[k]]$errVariance

    source("extractFunc.R")
    allocations = getAllocations(mcmcfit, mcmcIterNum = k)
    randComp = getRandomComp(mcmcfit, mcmcIterNum = k)
    
    Dthetabar = 0
    for(i in 1:nsubjects){
      startIndex = cumsumHb[i]+1
      endIndex = cumsumHb[i]+numHb[i]
      randomPartMean = randComp[i,1]*INTERCEPT_SCALE + randComp[i,2]*ds$donationLast2Years[startIndex:endIndex]

      fixedPartMean = betaAge*ds$Age[startIndex:endIndex] + betaSeason*(as.numeric(ds$Season[startIndex:endIndex])-1) +
        betaDonate*ds$Donate[startIndex:endIndex] + betaTSPD*ds$TSPD[startIndex:endIndex] +
        betaDonateLast2TSPD*ds$donateLast2TSPD[startIndex:endIndex] +
        betaDonateLast2Donate*ds$donateLast2Donate[startIndex:endIndex] +
        betaDonateLast2Square*ds$donateLast2Square[startIndex:endIndex]

       Dthetabar = Dthetabar + log(eta[allocations[i]]) +
         dmvnorm(x = ds$Hb[startIndex:endIndex],
                 mean = randomPartMean + fixedPartMean,
                 sigma = diag(numHb[i])*errVariance, log = T) +
         dmvnorm(x = randComp[i,],
                 mean = randmu[[allocations[i]]],
                 sigma = randSigma[[allocations[i]]], log = T)
    }
    Dthetabar
  }
  
  totalDthetabar = totalDthetabar[!(totalDthetabar %in% c(Inf,-Inf, NaN))]
  
  getPDandDIC(calculateCompleteMeanPostDeviance(mcmcfit), -2*mean(totalDthetabar))
}

calculate_E_theta_givenY_Z = function(mcmcfit){
  
  ncomponents = attributes(mcmcfit)$ncomponents
  nsubjects = attributes(mcmcfit)$nsubjects
  
  ds = ds
  cumsumHb = cumsumHb
  numHb = numHb
  INTERCEPT_SCALE=INTERCEPT_SCALE
  
  theta_givenY_Z=foreach(k=1:nrow(mcmcfit[[1]])) %dopar% {
    
    source("extractFunc.R")
    allocations = getAllocations(mcmcfit, mcmcIterNum = k)
    randComp = getRandomComp(mcmcfit, mcmcIterNum = k)
    
    #Eta
    freqPerGroup = sapply(1:ncomponents, function(j){sum(allocations==j)})
    eta=(1+freqPerGroup)/(ncomponents*1 + nsubjects)
    
    #randmu and #randSigma
    randmu = rep(list(c(0,0)), ncomponents)
    randSigma=rep(list(diag(2)),ncomponents)
    
    randCompPerGroup = lapply(1:ncomponents, function(x){matrix(randComp[allocations==x,], ncol=2)})
    
    for(j in 1:ncomponents){
      if(freqPerGroup[j]>0){
        randmu[[j]] = apply(randCompPerGroup[[j]], MARGIN=2, FUN=mean)
        randSigma[[j]] = var(randCompPerGroup[[j]])
      }
    }
    
    dsHb = ds$Hb
    dsAge = ds$Age
    dsDonate = as.numeric(ds$Donate)
    dsSeason = as.numeric(ds$Season)-1
    dsTSPD = ds$TSPD
    dsDonateLast2TSPD = ds$donateLast2TSPD
    dsDonateLast2Donate =  ds$donateLast2Donate
    dsDonateLast2Square = ds$donateLast2Square
    
    #Beta for regression
    for(i in 1:nsubjects){
      for(j in 1:numHb[i]){
        index = cumsumHb[i]+j
        dsHb[index] = dsHb[index] - randComp[i,1]*INTERCEPT_SCALE - randComp[i,2]*ds$donationLast2Years[index]
      }
    }
    
    reg=lm(dsHb~dsAge + dsDonate + dsSeason + dsTSPD + dsDonateLast2TSPD + dsDonateLast2Donate + dsDonateLast2Square + 0)
    errVariance = summary(reg)$sigma^2 * ((sum(numHb)-7)/(sum(numHb)-7-2))
    
    list("eta"=eta, "randmu"=randmu, "randSigma"=randSigma, 
         "betaAge"=reg$coefficients["dsAge"], "betaDonate"=reg$coefficients["dsDonate"], 
         "betaSeason"=reg$coefficients["dsSeason"], "betaTSPD" = reg$coefficients["dsTSPD"],
         "betaDonateLast2TSPD" = reg$coefficients["dsDonateLast2TSPD"],
         "betaDonateLast2Donate" = reg$coefficients["dsDonateLast2Donate"],
         "betaDonateLast2Square" = reg$coefficients["dsDonateLast2Square"],
         "errVariance"=errVariance)
  }
  
  return(theta_givenY_Z)
}

calculateCompleteMeanPostDeviance=function(mcmcfit){
  
  ds = ds
  cumsumHb = cumsumHb
  numHb = numHb
  INTERCEPT_SCALE = INTERCEPT_SCALE
  
  nsubjects = attributes(mcmcfit)$nsubjects
  
  totalDtheta=foreach(k=1:nrow(mcmcfit[[1]]), .combine = 'c', .packages = c('mvtnorm')) %dopar%{
    source("extractFunc.R")
    
    eta = getEta(mcmcfit, mcmcIterNum = k)
    
    betaAge = mcmcfit[[1]][k,"betaAge"]
    betaDonate = mcmcfit[[1]][k,"betaDonate"]
    betaSeason = mcmcfit[[1]][k,"betaSeason"]
    betaTSPD = mcmcfit[[1]][k,"betaTSPD"]
    betaDonateLast2TSPD = mcmcfit[[1]][k,"betaDonateLast2TSPD"]
    betaDonateLast2Donate = mcmcfit[[1]][k,"betaDonateLast2Donate"]
    betaDonateLast2Square = mcmcfit[[1]][k,"betaDonateLast2Square"]
    
    randComp = getRandomComp(mcmcfit, mcmcIterNum = k)
    randmu = getRandMu(mcmcfit, mcmcIterNum = k)
    errVariance = 1/mcmcfit[[1]][k,"errPrecision"]
    randSigma = getRandSigma(mcmcfit, mcmcIterNum = k)
    allocations = getAllocations(mcmcfit, mcmcIterNum = k)
    
    Dtheta = 0
    for(i in 1:nsubjects){
      startIndex = cumsumHb[i]+1
      endIndex = cumsumHb[i]+numHb[i]
      
      fixedPartMean = betaAge*ds$Age[startIndex:endIndex] + betaSeason*(as.numeric(ds$Season[startIndex:endIndex])-1) + 
        betaDonate*ds$Donate[startIndex:endIndex] + betaTSPD*ds$TSPD[startIndex:endIndex] + 
        betaDonateLast2TSPD*ds$donateLast2TSPD[startIndex:endIndex] +
        betaDonateLast2Donate*ds$donateLast2Donate[startIndex:endIndex] +
        betaDonateLast2Square*ds$donateLast2Square[startIndex:endIndex]
      
      randomPartMean = randComp[i,1]*INTERCEPT_SCALE + randComp[i,2]*ds$donationLast2Years[startIndex:endIndex]
      
      Dtheta = Dtheta + log(eta[allocations[i]]) + 
        dmvnorm(x = ds$Hb[startIndex:endIndex], 
                mean = randomPartMean + fixedPartMean, 
                sigma = diag(numHb[i])*errVariance, log = T) +
        dmvnorm(x = randComp[i,], 
                mean = randmu[[allocations[i]]], 
                sigma = randSigma[[allocations[i]]], log = T)
    }
    -2*Dtheta
  }
  mean(totalDtheta)
}