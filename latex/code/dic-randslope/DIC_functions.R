library(mvtnorm)
 
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
  
  ncomponents = attributes(mcmcfit)$ncomponents
  nsubjects = attributes(mcmcfit)$nsubjects

  totalDtheta=foreach(k=1:nrow(mcmcfit[[1]]), .combine = c, .packages = c('mvtnorm')) %dopar%{
  
    source("../common/extractFuncRandSlopeConditional.R")

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

calculateDIC1 = function(mcmcfit){
  getPDandDIC(calculateObsMeanPostDeviance(mcmcfit), calculateObsDThetaBar(mcmcfit, "mean"))
}

calculateDIC2 = function(mcmcfit){
  getPDandDIC(calculateObsMeanPostDeviance(mcmcfit), calculateObsDThetaBar(mcmcfit, "modeFunc"))
}

calculateDIC3 = function(mcmcfit){
  
  DThetabar_functionalapprox = 0
  
  dswide = dswide
  dswide_y = dswide_y
  time = time
  
  ncomponents = attributes(mcmcfit)$ncomponents
  nsubjects = attributes(mcmcfit)$nsubjects

  obsDev_Subjects=foreach(k=1:nrow(mcmcfit[[1]]), .combine = 'rbind', .packages = c('mvtnorm')) %dopar%{
    
    source("../common/extractFuncRandSlopeConditional.R")
    
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

calculateDIC7 = function(mcmcfit, summaryFuncName="modeFunc"){
  meanPostDeviance = 0
  
  dswide = dswide
  dswide_y = dswide_y
  time = time
  
  ncomponents = attributes(mcmcfit)$ncomponents
  nsubjects = attributes(mcmcfit)$nsubjects
  nrep = attributes(mcmcfit)$nrep
  
  totalDtheta=foreach(k=1:nrow(mcmcfit[[1]]), .combine = 'c', .packages = c('mvtnorm')) %dopar%{
    source("../common/extractFuncRandSlopeConditional.R")
    
    betaBy = mcmcfit[[1]][k,"betaBy"]
    betaGender = mcmcfit[[1]][k,"betaGender"]
    randComp = getRandomComp(mcmcfit, mcmcIterNum = k)
    errVariance = 1/mcmcfit[[1]][k,"errPrecision"]
    
    Dtheta = 0
    for(i in 1:nsubjects){
      randomPartMean = randComp[i,1] + randComp[i,2]*time
      fixedPartMean = betaBy * (as.numeric(dswide[i, "by"])-1) + 
        betaGender * (as.numeric(dswide[i, "gender"])-1)
      Dtheta = Dtheta + dmvnorm(x = dswide_y[i,], 
                                     mean = randomPartMean + fixedPartMean, 
                                     sigma = diag(nrep)*errVariance, log = T)
    }
    Dtheta
  }
  meanPostDeviance = -2*mean(totalDtheta)
  
  #Calcuing D(thetabar)
  summaryFunc = get(summaryFuncName)
  
  randComp = getRandomComp(mcmcfit, summaryFuncName)
  betaBy = summaryFunc(mcmcfit[[1]][,"betaBy"])
  betaGender = summaryFunc(mcmcfit[[1]][,"betaGender"])
  errVariance = summaryFunc(1/mcmcfit[[1]][,"errPrecision"])
  
  Dthetabar = 0
  for(i in 1:nsubjects){
    randomPartMean = randComp[i,1] + randComp[i,2]*time
    fixedPartMean = betaBy * (as.numeric(dswide[i, "by"])-1) + 
      betaGender * (as.numeric(dswide[i, "gender"])-1)
    Dthetabar = Dthetabar + dmvnorm(x = dswide_y[i,], 
                              mean = randomPartMean + fixedPartMean, 
                              sigma = diag(nrep)*errVariance, log = T)
  }
  
  getPDandDIC(meanPostDeviance, -2*Dthetabar)
}

calculateDIC5 = function(mcmcfit){
  summaryFuncName = "modeFunc"
  
  #Calcuing D(thetabar)
  summaryFunc = get(summaryFuncName)
  ncomponents = attributes(mcmcfit)$ncomponents
  nsubjects = attributes(mcmcfit)$nsubjects
  nrep = attributes(mcmcfit)$nrep
  
  eta = getEta(mcmcfit, summaryFuncName)
  randComp = getRandomComp(mcmcfit, summaryFuncName)
  randmu = getRandMu(mcmcfit, summaryFuncName)
  randSigma = getRandSigma(mcmcfit, summaryFuncName)
  betaBy = summaryFunc(mcmcfit[[1]][,"betaBy"])
  betaGender = summaryFunc(mcmcfit[[1]][,"betaGender"])
  errVariance = summaryFunc(1/mcmcfit[[1]][,"errPrecision"])
  allocations = getAllocations(mcmcfit, summaryFuncName)
  
  Dthetabar = 0
  for(i in 1:nsubjects){
    randomPartMean = randComp[i,1] + randComp[i,2]*time
    fixedPartMean = betaBy * (as.numeric(dswide[i, "by"])-1) + 
      betaGender * (as.numeric(dswide[i, "gender"])-1)
    Dthetabar = Dthetabar + log(eta[allocations[i]])+
      dmvnorm(x = dswide_y[i,], 
              mean = randomPartMean + fixedPartMean, 
              sigma = diag(nrep)*errVariance, log = T)+
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
  nrep = attributes(mcmcfit)$nrep
     
  #Calcuing D(thetabar)
  totalDthetabar = foreach(k=1:nrow(mcmcfit[[1]]), .combine = 'c', .packages = c('mvtnorm')) %dopar%{
    eta = E_theta_givenY_Z[[k]]$eta
    randmu = E_theta_givenY_Z[[k]]$randmu
    randSigma = E_theta_givenY_Z[[k]]$randSigma
    betaBy = E_theta_givenY_Z[[k]]$betaBy
    betaGender = E_theta_givenY_Z[[k]]$betaGender
    errVariance = E_theta_givenY_Z[[k]]$errVariance

    allocations = getAllocations(mcmcfit, mcmcIterNum = k)
    randComp = getRandomComp(mcmcfit, mcmcIterNum = k)
    
    Dthetabar = 0
    for(i in 1:nsubjects){
      randomPartMean = randComp[i,1] + randComp[i,2]*time
      fixedPartMean = betaBy * (as.numeric(dswide[i, "by"])-1) + 
        betaGender * (as.numeric(dswide[i, "gender"])-1)
      Dthetabar = Dthetabar + log(eta[allocations[i]])+
        dmvnorm(x = dswide_y[i,], 
                mean = randomPartMean + fixedPartMean, 
                sigma = diag(nrep)*errVariance, log = T)+
        dmvnorm(x = randComp[i,], 
                mean = randmu[[allocations[i]]], 
                sigma = randSigma[[allocations[i]]], log = T)
    }
    Dthetabar
  }
  
  getPDandDIC(calculateCompleteMeanPostDeviance(mcmcfit), -2*mean(totalDthetabar))
}

calculate_E_theta_givenY_Z = function(mcmcfit){
  
  ncomponents = attributes(mcmcfit)$ncomponents
  nsubjects = attributes(mcmcfit)$nsubjects
  nrep = attributes(mcmcfit)$nrep
  dsweight = ds$weight
  dstime = ds$time
  dsgender = as.numeric(ds$gender)-1
  dsby = as.numeric(ds$by)-1
  
  theta_givenY_Z=foreach(k=1:nrow(mcmcfit[[1]])) %dopar% {
    
    source("../common/extractFuncRandSlopeConditional.R")
    allocations = getAllocations(mcmcfit, mcmcIterNum = k)
    randComp = getRandomComp(mcmcfit, mcmcIterNum = k)
    
    #Eta
    freqPerGroup = rep(0, ncomponents)
    tempFreq = table(allocations)
    for(j in dimnames(tempFreq)[[1]]){
      freqPerGroup[as.numeric(j)] = tempFreq[j]
    }
    eta=(1+freqPerGroup)/(ncomponents + nsubjects)
    
    #randmu and #randSigma
    wishartPriorDf = 3
    randmu = rep(list(c(0,0)), ncomponents)
    randSigma=rep(list(diag(2)),ncomponents)
    randCompPerGroup = rep(list(matrix(nrow=0, ncol=2)), ncomponents)
    for(j in 1:nsubjects){
      randCompPerGroup[[allocations[j]]] = rbind(randCompPerGroup[[allocations[j]]], randComp[j,])
    }
    for(j in 1:ncomponents){
      if(freqPerGroup[j]>0){
        randmu[[j]] = apply(randCompPerGroup[[j]], MARGIN=2, FUN=mean)
        
        varIntercept = (sum((randCompPerGroup[[j]][,1]-randmu[[j]][1])^2)+0.001)/(freqPerGroup[j]+0.001-2)
        varSlope = (sum((randCompPerGroup[[j]][,2]-randmu[[j]][2])^2)+0.001)/(freqPerGroup[j]+0.001-2)
        covIntSlope = cov(randCompPerGroup[[j]][,1], randCompPerGroup[[j]][,2])

        randSigma[[j]] = matrix(c(varIntercept, covIntSlope, covIntSlope, varSlope), nrow=2, ncol=2)
      }
    }
    
    #Beta for regression
    for(i in 1:nsubjects){
      for(j in 1:nrep){
        index = (i-1)*nrep + j
        dsweight[index] = dsweight[index] - randComp[i,1] - randComp[i,2]*dstime[index]
      }
    }
    
    reg=lm(dsweight~dsgender + dsby + 0)
    errVariance = summary(reg)$sigma^2
    
    list("eta"=eta, "randmu"=randmu, "randSigma"=randSigma, 
         "betaBy"=reg$coefficients["dsby"], 
         "betaGender"=reg$coefficients["dsgender"], "errVariance"=errVariance)
  }
  
  return(theta_givenY_Z)
}

calculateCompleteMeanPostDeviance=function(mcmcfit){
  
  dswide = dswide
  dswide_y = dswide_y
  time = time
  
  nsubjects = attributes(mcmcfit)$nsubjects
  nrep = attributes(mcmcfit)$nrep
  
  totalDtheta=foreach(k=1:nrow(mcmcfit[[1]]), .combine = 'c', .packages = c('mvtnorm')) %dopar%{
    source("../common/extractFuncRandSlopeConditional.R")
    
    eta = getEta(mcmcfit, mcmcIterNum = k)
    betaBy = mcmcfit[[1]][k,"betaBy"]
    betaGender = mcmcfit[[1]][k,"betaGender"]
    randComp = getRandomComp(mcmcfit, mcmcIterNum = k)
    randmu = getRandMu(mcmcfit, mcmcIterNum = k)
    errVariance = 1/mcmcfit[[1]][k,"errPrecision"]
    randSigma = getRandSigma(mcmcfit, mcmcIterNum = k)
    allocations = getAllocations(mcmcfit, mcmcIterNum = k)
    
    Dtheta = 0
    for(i in 1:nsubjects){
      randomPartMean = randComp[i,1] + randComp[i,2]*time
      fixedPartMean = betaBy * (as.numeric(dswide[i, "by"])-1) + 
        betaGender * (as.numeric(dswide[i, "gender"])-1)
      
      Dtheta = Dtheta + log(eta[allocations[i]]) + 
        dmvnorm(x = dswide_y[i,], 
                mean = randomPartMean + fixedPartMean, 
                sigma = diag(nrep)*errVariance, log = T) +
        dmvnorm(x = randComp[i,], 
                mean = randmu[[allocations[i]]], 
                sigma = randSigma[[allocations[i]]], log = T)
    }
    -2*Dtheta
  }
  mean(totalDtheta)
}