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
  betaAge = summaryFunc(mcmcfit[[1]][,"betaAge"])
  
  Dthetabar = 0
  for(i in 1:nsubjects){
    temp = 0
    for(j in 1:ncomponents){
      randomPartMean = randmuSummary[[j]][1] + randmuSummary[[j]][2]*time
      fixedPartMean = betaBy * (as.numeric(dswide[i, "by"])-1) + 
        betaGender * (as.numeric(dswide[i, "gender"])-1) + betaAge * dswide[i,"age"]
      temp = temp + etaSummary[j] * dmvnorm(x = dswide_y[i,], 
                                            mean = randomPartMean + fixedPartMean, 
                                            sigma = sigmaSummary[[j]], log = F)
    }
    Dthetabar = Dthetabar + log(temp)
  }
  Dthetabar = -2*Dthetabar
}

calculateObsTotalDTheta = function(mcmcfit){
  dswide = dswide
  dswide_y = dswide_y
  time = attributes(mcmcfit)$time
  
  ncomponents = attributes(mcmcfit)$ncomponents
  nsubjects = attributes(mcmcfit)$nsubjects
  
  totalDtheta=foreach(k=1:nrow(mcmcfit[[1]]), .combine = c, .packages = c('mvtnorm')) %dopar%{
   #for(k in 1:nrow(mcmcfit[1])){ 
    
    source("../common/extractFuncRandSlopeConditional.R")
    eta = getEta(mcmcfit, mcmcIterNum = k)
    sigma = getSigma(mcmcfit, mcmcIterNum = k)
    randmu = getRandMu(mcmcfit, mcmcIterNum = k)
    betaBy = mcmcfit[[1]][k,"betaBy"]
    betaGender = mcmcfit[[1]][k,"betaGender"]
    betaAge = mcmcfit[[1]][k,"betaAge"]
    
    Dtheta = 0
    for(i in 1:nsubjects){
      temp = 0
      for(j in 1:ncomponents){
        randomPartMean = randmu[[j]][1] + randmu[[j]][2]*time
        fixedPartMean = betaBy * (as.numeric(dswide[i, "by"])-1) + 
          betaGender * (as.numeric(dswide[i, "gender"])-1) + betaAge * dswide[i, "age"]
        temp = temp + eta[j] * dmvnorm(x = dswide_y[i,], 
                                       mean = randomPartMean + fixedPartMean, 
                                       sigma = sigma[[j]], log = F)
      }
      Dtheta = Dtheta + log(temp)
    }
    -2*Dtheta
  }
  totalDtheta
}

calculateDIC1 = function(mcmcfit){
  getPDandDIC(mean(calculateObsTotalDTheta(mcmcfit)), calculateObsDThetaBar(mcmcfit, "mean"))
}

calculateDIC2 = function(mcmcfit){
  getPDandDIC(mean(calculateObsTotalDTheta(mcmcfit)), calculateObsDThetaBar(mcmcfit, "modeFunc"))
}

calculateDIC3 = function(mcmcfit){
  
  DThetabar_functionalapprox = 0
  
  dswide = dswide
  dswide_y = dswide_y
  time = attributes(mcmcfit)$time
  
  ncomponents = attributes(mcmcfit)$ncomponents
  nsubjects = attributes(mcmcfit)$nsubjects

  obsDev_Subjects=foreach(k=1:nrow(mcmcfit[[1]]), .combine = 'rbind', .packages = c('mvtnorm')) %dopar%{
    
    source("../common/extractFuncRandSlopeConditional.R")
    
    eta = getEta(mcmcfit, mcmcIterNum = k)
    sigma = getSigma(mcmcfit, mcmcIterNum = k)
    randmu = getRandMu(mcmcfit, mcmcIterNum = k)
    betaBy = mcmcfit[[1]][k,"betaBy"]
    betaGender = mcmcfit[[1]][k,"betaGender"]
    betaAge = mcmcfit[[1]][k,"betaAge"]
    
    Dtheta = rep(0, nsubjects)
    for(i in 1:nsubjects){
      for(j in 1:ncomponents){
        randomPartMean = randmu[[j]][1] + randmu[[j]][2]*time
        fixedPartMean = betaBy * (as.numeric(dswide[i, "by"])-1) + 
          betaGender * (as.numeric(dswide[i, "gender"])-1) + betaAge * dswide[i, "age"]
        Dtheta[i] = Dtheta[i] + eta[j] * dmvnorm(x = dswide_y[i,], 
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
  
  dswide = dswide
  dswide_y = dswide_y
  time = attributes(mcmcfit)$time
  
  ncomponents = attributes(mcmcfit)$ncomponents
  nsubjects = attributes(mcmcfit)$nsubjects
  nrep = attributes(mcmcfit)$nrep
  
  totalDtheta=foreach(k=1:nrow(mcmcfit[[1]]), .combine = 'c', .packages = c('mvtnorm')) %dopar%{
    source("../common/extractFuncRandSlopeConditional.R")
    
    betaBy = mcmcfit[[1]][k,"betaBy"]
    betaGender = mcmcfit[[1]][k,"betaGender"]
    betaAge = mcmcfit[[1]][k,"betaAge"]
    randComp = getRandomComp(mcmcfit, mcmcIterNum = k)
    errVariance = 1/mcmcfit[[1]][k,"errPrecision"]
    
    Dtheta = 0
    for(i in 1:nsubjects){
      randomPartMean = randComp[i,1] + randComp[i,2]*time
      fixedPartMean = betaBy * (as.numeric(dswide[i, "by"])-1) + 
        betaGender * (as.numeric(dswide[i, "gender"])-1) + betaAge * dswide[i, "age"]
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
  betaAge = summaryFunc(mcmcfit[[1]][,"betaAge"])
  errVariance = summaryFunc(1/mcmcfit[[1]][,"errPrecision"])
  
  Dthetabar = 0
  for(i in 1:nsubjects){
    randomPartMean = randComp[i,1] + randComp[i,2]*time
    fixedPartMean = betaBy * (as.numeric(dswide[i, "by"])-1) + 
      betaGender * (as.numeric(dswide[i, "gender"])-1) + betaAge * dswide[i, "age"]
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
  betaAge = summaryFunc(mcmcfit[[1]][,"betaAge"])
  errVariance = summaryFunc(1/mcmcfit[[1]][,"errPrecision"])
  allocations = getAllocations(mcmcfit, summaryFuncName)
  
  Dthetabar = 0
  for(i in 1:nsubjects){
    randomPartMean = randComp[i,1] + randComp[i,2]*time
    fixedPartMean = betaBy * (as.numeric(dswide[i, "by"])-1) + 
      betaGender * (as.numeric(dswide[i, "gender"])-1) + betaAge * dswide[i, "age"]
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
     
  dswide = dswide
  dswide_y = dswide_y
  time = attributes(mcmcfit)$time
  
  #Calcuing D(thetabar)
  totalDthetabar = foreach(k=1:nrow(mcmcfit[[1]]), .combine = 'c', .packages = c('mvtnorm')) %dopar%{
    eta = E_theta_givenY_Z[[k]]$eta
    randmu = E_theta_givenY_Z[[k]]$randmu
    randSigma = E_theta_givenY_Z[[k]]$randSigma
    betaBy = E_theta_givenY_Z[[k]]$betaBy
    betaGender = E_theta_givenY_Z[[k]]$betaGender
    betaAge = E_theta_givenY_Z[[k]]$betaAge
    errVariance = E_theta_givenY_Z[[k]]$errVariance

    source("../common/extractFuncRandSlopeConditional.R")
    allocations = getAllocations(mcmcfit, mcmcIterNum = k)
    randComp = getRandomComp(mcmcfit, mcmcIterNum = k)
    
    Dthetabar = 0
    for(i in 1:nsubjects){
      randomPartMean = randComp[i,1] + randComp[i,2]*time
       fixedPartMean = betaBy * (as.numeric(dswide[i, "by"])-1) + 
         betaGender * (as.numeric(dswide[i, "gender"])-1) + betaAge*dswide[i, "age"]
       Dthetabar = Dthetabar + log(eta[allocations[i]])+
         dmvnorm(x = dswide_y[i,], 
                 mean = randomPartMean + fixedPartMean, 
                 sigma = diag(nrep)*errVariance, log = T)
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
  ds=ds
  
  theta_givenY_Z=foreach(k=1:nrow(mcmcfit[[1]])) %dopar% {
    
    source("../common/extractFuncRandSlopeConditional.R")
    allocations = getAllocations(mcmcfit, mcmcIterNum = k)
    randComp = getRandomComp(mcmcfit, mcmcIterNum = k)
    
    #Eta
    freqPerGroup = sapply(1:ncomponents, function(j){sum(allocations==j)})
    eta=(3+freqPerGroup)/(ncomponents*3 + nsubjects)
    
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
    
    dsweight = ds$weight
    dstime = ds$time
    dsgender = as.numeric(ds$gender)-1
    dsby = as.numeric(ds$by)-1
    dsage = ds$age
    
    #Beta for regression
    for(i in 1:nsubjects){
      for(j in 1:nrep){
        index = (i-1)*nrep + j
        dsweight[index] = dsweight[index] - randComp[i,1] - randComp[i,2]*dstime[index]
      }
    }
    
    reg=lm(dsweight~dsgender + dsby + dsage + 0)
    errVariance = summary(reg)$sigma^2 * ((nsubjects*nrep-3)/(nsubjects*nrep-3-2))
    
    list("eta"=eta, "randmu"=randmu, "randSigma"=randSigma, 
         "betaBy"=reg$coefficients["dsby"], 
         "betaGender"=reg$coefficients["dsgender"], "betaAge"=reg$coefficients["dsage"], "errVariance"=errVariance)
  }
  
  return(theta_givenY_Z)
}

calculateCompleteMeanPostDeviance=function(mcmcfit){
  
  dswide = dswide
  dswide_y = dswide_y
  time = attributes(mcmcfit)$time
  
  nsubjects = attributes(mcmcfit)$nsubjects
  nrep = attributes(mcmcfit)$nrep
  
  totalDtheta=foreach(k=1:nrow(mcmcfit[[1]]), .combine = 'c', .packages = c('mvtnorm')) %dopar%{
    source("../common/extractFuncRandSlopeConditional.R")
    
    eta = getEta(mcmcfit, mcmcIterNum = k)
    betaBy = mcmcfit[[1]][k,"betaBy"]
    betaGender = mcmcfit[[1]][k,"betaGender"]
    betaAge = mcmcfit[[1]][k,"betaAge"]
    
    randComp = getRandomComp(mcmcfit, mcmcIterNum = k)
    randmu = getRandMu(mcmcfit, mcmcIterNum = k)
    errVariance = 1/mcmcfit[[1]][k,"errPrecision"]
    randSigma = getRandSigma(mcmcfit, mcmcIterNum = k)
    allocations = getAllocations(mcmcfit, mcmcIterNum = k)
    
    Dtheta = 0
    for(i in 1:nsubjects){
      randomPartMean = randComp[i,1] + randComp[i,2]*time
      fixedPartMean = betaBy * (as.numeric(dswide[i, "by"])-1) + 
        betaGender * (as.numeric(dswide[i, "gender"])-1) + betaAge * dswide[i, "age"]
      
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