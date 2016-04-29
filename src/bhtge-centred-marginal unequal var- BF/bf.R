ALMOST_ZERO = 0.00000000000001

modeFunc=function(data){
  d = density(data, n=1e4)
  return(d$x[which.max(d$y)])
}

getRandomIntercept=function(mcmcfit, mcmcIterNum){
  randomIntercept = numeric()
  for(i in 1:nsubjects){
    randomIntercept[i] = mcmcfit[[1]][mcmcIterNum, paste("randomIntercept[",i,"]", sep="")]
  }
  return(randomIntercept)
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

getRandMu=function(mcmcfit, mcmcIterNum = NA, summaryFuncName="mean"){
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

getRandPrecision=function(mcmcfit, mcmcIterNum = NA, summaryFuncName="mean"){
  summaryFunc = get(summaryFuncName)
  randPrecision=numeric()
  
  for(i in 1:ncomponents){
    if(!is.numeric(mcmcIterNum)){
      randPrecision[i] = summaryFunc(mcmcfit[[1]][,paste("randPrecision[",i,"]", sep="")])
    }else{
      randPrecision[i] = mcmcfit[[1]][mcmcIterNum, paste("randPrecision[",i,"]", sep="")]
    }
  }
  return(randPrecision + ALMOST_ZERO)
}

getSigma=function(mcmcfit, mcmcIterNum = NA, componentNum=1){
  randVar = 1/(mcmcfit[[1]][mcmcIterNum, paste("randPrecision[",componentNum,"]", sep="")] + ALMOST_ZERO)
  errVar = 1/mcmcfit[[1]][mcmcIterNum, "errPrecision"]
  
  sigma = diag(nrep) * errVar
  for(i in 1:nrep){
    for(j in 1:nrep){
      sigma[i,j] = sigma[i,j] + randVar
    }
  }

  return(sigma)
}

getAllSigma = function(mcmcfit, mcmcIterNum=NA, ncomp){
  allSigma = list(ncomp)
  for(i in 1:ncomp){
    allSigma[[i]] = getSigma(mcmcfit, mcmcIterNum, i)
  }
  return(allSigma)
}

getXBeta=function(mcmcfit, mcmcIterNum = NA){
  XBeta=matrix(data=NA, nrow=nsubjects, ncol = nrep)
  
  for(i in 1:nsubjects){
    for(j in 1:nrep){
      XBeta[i,j] = mcmcfit[[1]][mcmcIterNum, paste("XBeta[",(i-1)*nrep + j,"]", sep="")]
    }
  }

  return(XBeta)
}

getAllocations=function(mcmcfit, mcmcIterNum){
  allocation = numeric(nsubjects)
  for(j in 1:nsubjects){
    allocation[j] = mcmcfit[[1]][mcmcIterNum, paste("S[",j,"]", sep="")]
  }
  return(allocation)
}

#Step 1: logL
logL = function(datawide, eta, sigma, bGender, bBy, bTime, randmu){
  weights_wide=as.matrix(dswide[,c(-1,-2,-3)])
  gender_wide = as.numeric(datawide$gender)-1
  birthyear_wide = as.numeric(datawide$by)-1
  nrep = ncol(datawide)-3
  
  logLikelihood = 0
  for(i in 1:nsubjects){
    temp = 0
    for(j in 1:ncomponents){
      tempmean = bTime*(1:nrep) + bGender*gender_wide[i]+bBy*birthyear_wide[i] + randmu[j]
      temp = temp + eta[j] * dmvnorm(weights_wide[i,], mean=tempmean, sigma=sigma[[j]],log = F)
    }
    logLikelihood = logLikelihood + log(temp)
  }
  return(logLikelihood)
}

#Step 2: log prior p(theta*)
log_prior=function(eta, randmu, betaGender, betaBy, betaTime, errPrecision, randPrecision) {
  dirichParm=getDirichParam(ncomponents)
  lgamma(sum(dirichParm)) - sum(lgamma(dirichParm)) + sum((dirichParm-1)*log(eta)) + 
  sum(dnorm(c(randmu, betaGender, betaBy, betaTime), mean = betaMu, sd = sqrt(1/betaTau), log = T)) +
  sum(dgamma(c(errPrecision, randPrecision), shape = gammaShapeRate, rate = gammaShapeRate, log = T))
}

maxValueIter = which.max(foreach(i=1:mcmcLen, .combine = c, .packages = c('mvtnorm')) %dopar%{
  logL(dswide, getEta(mcmcfit, mcmcIterNum = i), getAllSigma(mcmcfit, mcmcIterNum = i,  ncomp = ncomponents),
         mcmcfit[[1]][i, "betaGender"],
         mcmcfit[[1]][i, "betaBy"],
         mcmcfit[[1]][i, "betaTime"], getRandMu(mcmcfit, mcmcIterNum = i))
  })

randmu_max = getRandMu(mcmcfit, mcmcIterNum = maxValueIter)
betaGender_max = mcmcfit[[1]][maxValueIter, "betaGender"]
betaBy_max = mcmcfit[[1]][maxValueIter, "betaBy"]
betaTime_max = mcmcfit[[1]][maxValueIter, "betaTime"]
errPrecision_max = mcmcfit[[1]][maxValueIter, "errPrecision"]
randPrecision_max = getRandPrecision(mcmcfit, mcmcIterNum = maxValueIter)
sigma_max = getAllSigma(mcmcfit, mcmcIterNum = maxValueIter, ncomp = ncomponents)
allocations_max = getAllocations(mcmcfit, mcmcIterNum = maxValueIter)
randIntercept_max = getRandomIntercept(mcmcfit, maxValueIter)
eta_max = getEta(mcmcfit, mcmcIterNum = maxValueIter)

chib = logL(dswide, eta_max, sigma_max, betaGender_max, betaBy_max, betaTime_max, 
            randmu_max) + log_prior(eta_max,randmu_max, betaGender_max, 
            betaBy_max, betaTime_max, errPrecision_max, randPrecision_max)

###### Conditional posterior of randprecision #######
#p(randPrecision|data)
mcmcfit_randPrecision = mcmcfit
chib = chib - log(mean(foreach(i=1:mcmcLen, .combine = c) %dopar%{
  
  randComp = getRandomIntercept(mcmcfit_randPrecision, i)
  allocations = getAllocations(mcmcfit_randPrecision,i)
  randmu = getRandMu(mcmcfit_randPrecision, i)

    meanCentredRandEff=numeric()
  for(j in 1:ncomponents){
    meanCentredRandEff = c(meanCentredRandEff, randComp[allocations==j]-randmu[j])
  }
  dgamma(randPrecision_max, shape = gammaShapeRate+nsubjects/2, rate=gammaShapeRate+sum(meanCentredRandEff^2)/2)
}))

######## Conditional posterior of randmu #########
#p(randMu|randPrecision,x)
mcmcfit_randmu = fitModel_randmu(niter, nthin, nburnin, jagsmodel = model_randmu, nchains = numchains, ncomponents,
                                 randPrecision_max)$mcmcfit

chib = chib - log(mean(foreach(i=1:mcmcLen, .combine = c) %dopar%{
  
  randComp = getRandomIntercept(mcmcfit_randmu, i)
  allocations = getAllocations(mcmcfit_randmu,i)
  
  n <- sapply(1:ncomponents, function(j){sum(allocations==j)})
  
  prob = 1
  for(j in 1:ncomponents){
    newMean = (sum(randComp[allocations==j])*randPrecision_max[j] + 
                 betaTau*betaMu)/(betaTau+n[j]*randPrecision_max[j])
    prob = prob * dnorm(randmu_max[j], mean=newMean, sd = 1/sqrt(betaTau + n[j]*randPrecision_max[j]))
  }
  return(prob)
}))

######## Conditional posterior of errPrecision #########
#p(errPrecision|randPrecision, randmu, x)
mcmcfit_errPrecision = fitModel_errPrecision(niter, nthin, nburnin, jagsmodel = model_errPrecision,
                    nchains = numchains, ncomponents, randmu_max, randPrecision_max)$mcmcfit

chib = chib - log(mean(foreach(i=1:mcmcLen, .combine = c) %dopar%{

    randComp = getRandomIntercept(mcmcfit_errPrecision, i)
    allocations = getAllocations(mcmcfit_errPrecision,i)
    bGender = mcmcfit_errPrecision[[1]][i, "betaGender"]
    bBy = mcmcfit_errPrecision[[1]][i, "betaBy"]
    bTime = mcmcfit_errPrecision[[1]][i, "betaTime"]
  
    n <- sapply(1:ncomponents, function(j){sum(allocations==j)})
    temp = ds$weight - bGender*(as.numeric(ds$gender)-1) - bBy*(as.numeric(ds$by)-1) - bTime*ds$time
    for(j in 1:nsubjects){
      startIndex = (j-1)*nrep + 1
      endIndex = j*nrep
      temp[startIndex:endIndex] = temp[startIndex:endIndex] - randComp[j]
    }
    
    dgamma(errPrecision_max, shape = gammaShapeRate+nsubjects*nrep/2, 
           rate=gammaShapeRate+sum(temp^2)/2)
}))

######## Conditional posterior of beta #########
#p(beta|errPrecision, randPrecision, randmu, x)
mcmcfit_beta = fitModel_beta(niter, nthin, nburnin, jagsmodel = model_beta, 
                         nchains = numchains, ncomponents, randmu_max, 
                         randPrecision_max, errPrecision_max)$mcmcfit

chib = chib - log(mean(foreach(i=1:mcmcLen, .combine = c) %dopar%{
  
  randComp = getRandomIntercept(mcmcfit_beta, i)
  allocations = getAllocations(mcmcfit_beta,i)
  
  n <- sapply(1:ncomponents, function(j){sum(allocations==j)})

  temp = ds$weight
  for(j in 1:nsubjects){
    startIndex = (j-1)*nrep + 1
    endIndex = j*nrep
    temp[startIndex:endIndex] = temp[startIndex:endIndex] - randComp[j]
  }
  
  tempreg = lm(temp~I(as.numeric(ds$gender)-1)+I(as.numeric(ds$by)-1)+ds$time+0)
  designMatrix = model.matrix(tempreg)
  dmvnorm(c(betaGender_max, betaBy_max, betaTime_max), mean = tempreg$coefficients, 
          sigma = solve(t(designMatrix)%*%designMatrix)/errPrecision_max)
}))

######## Conditional posterior of eta #########
#p(eta|errPrecision, randPrecision, randmu, beta, x)
mcmcfit_eta = fitModel_eta(niter, nthin, nburnin, jagsmodel = model_eta, 
                                             nchains = numchains, ncomponents, randmu_max, 
                                             randPrecision_max, errPrecision_max, betaGender_max,
                                            betaBy_max, betaTime_max)$mcmcfit

chib = chib - log(mean(foreach(i=1:mcmcLen, .combine = c) %dopar%{
  
  allocations = getAllocations(mcmcfit_eta,i)
  dirichParm=getDirichParam(ncomponents)
  n <- sapply(1:ncomponents, function(j){sum(allocations==j)})
  
  exp(lgamma(sum(dirichParm)+nsubjects) - sum(lgamma(dirichParm+n)) + 
        sum((dirichParm+n-1)*log(eta_max)))
}))



