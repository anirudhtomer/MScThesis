#Based on Chib's approximation

source("DIC_functions.R")

#Step 1: logL
#calculated via the DIC function calculateObsTotalDTheta
#since the function I am calling gives -2*logL, I multiply it with -0.5
obslogL=calculateObsTotalDTheta(mcmcfit)*(-0.5)
maxValueIter = which.max(obslogL)

randmu_max = getRandMu(mcmcfit, mcmcIterNum = maxValueIter)
randIntercept_max = sapply(randmu_max, function(x){x[1]})
randSlope_max = sapply(randmu_max, function(x){x[2]})

betaGender_max = mcmcfit[[1]][maxValueIter, "betaGender"]
betaBy_max = mcmcfit[[1]][maxValueIter, "betaBy"]

errPrecision_max = mcmcfit[[1]][maxValueIter, "errPrecision"]

randPrecision_max = lapply(getRandSigma(mcmcfit, mcmcIterNum = maxValueIter), FUN = function(sigmaMatrix){solve(sigmaMatrix)})
precisionIntercept_max = sapply(randPrecision_max, FUN = function(mat){mat[1,1]})
precisionSlope_max = sapply(randPrecision_max, FUN = function(mat){mat[2,2]})

sigma_max = getSigma(mcmcfit, mcmcIterNum = maxValueIter)
eta_max = getEta(mcmcfit, mcmcIterNum = maxValueIter)

#Step 2: log prior p(theta*)
log_prior=lgamma(sum(dirichParm)) - sum(lgamma(dirichParm)) + sum((dirichParm-1)*log(eta_max)) + 
  sum(dnorm(c(randIntercept_max, randSlope_max, betaGender_max, betaBy_max), 
            mean = betaMu, sd = sqrt(1/betaTau), log = T)) +
  sum(dgamma(c(errPrecision_max, precisionIntercept_max, precisionSlope_max), 
             shape = gammaShapeRate, rate = gammaShapeRate, log = T))

chib = obslogL[maxValueIter] + log_prior

###### Conditional posterior of randprecision #######
#p(randPrecision|data)
mcmcfit_randPrecision = mcmcfit
chib = chib - log(mean(foreach(i=1:nrow(mcmcfit[[1]]), .combine = c) %dopar%{
  
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
  
  dgamma(randPrecision_max, shape = gammaShapeRate+nsubjects/2, rate=gammaShapeRate+sum(meanCentredRandEff^2)/2)
}))

######## Conditional posterior of randmu #########
#p(randMu|randPrecision,x)
mcmcfit_randmu = fitModel_randmu(niter, nthin, nburnin, jagsmodel = model_randmu, nchains = numchains, ncomponents,
                                 randPrecision_max)$mcmcfit

chib = chib - log(mean(foreach(i=1:nrow(mcmcfit_randmu[[1]]), .combine = c) %dopar%{
  ncomponents=attributes(mcmcfit_randmu)$ncomponents
  
  randComp = getRandomComp(mcmcfit_randmu, i)
  allocations = getAllocations(mcmcfit_randmu, i)
  
  nPerGroup <- sapply(1:ncomponents, function(j){sum(allocations==j)})
  
  prob = 1
  for(j in 1:ncomponents){
    newMean_Intercept = (sum(randComp[allocations==j])*randPrecision_max + betaTau*betaMu)/(betaTau+n[j]*randPrecision_max)
    prob = prob * dnorm(randmu_max[j], mean=newMean, sd = 1/sqrt(betaTau+n[j]*randPrecision_max))
  }
  return(prob)
}))

######## Conditional posterior of errPrecision #########
#p(errPrecision|randPrecision, randmu, x)
mcmcfit_errPrecision = fitModel_errPrecision(niter, nthin, nburnin, jagsmodel = model_errPrecision,
                    nchains = numchains, ncomponents, randmu_max, randPrecision_max)$mcmcfit

chib = chib - log(mean(foreach(i=1:nrow(mcmcfit_errPrecision[[1]]), .combine = c) %dopar%{

    randComp = getRandomComp(mcmcfit_errPrecision, i)
    allocations = getAllocations(mcmcfit_errPrecision,i)
    bGender = mcmcfit_errPrecision[[1]][i, "betaGender"]
    bBy = mcmcfit_errPrecision[[1]][i, "betaBy"]
    
    ncomponents = attributes(mcmcfit_errPrecision)$ncomponents
    nsubjects = attributes(mcmcfit_errPrecision)$nsubjects
    nrep = attributes(mcmcfit_errPrecision)$nrep
    time = attributes(mcmcfit_errPrecision)$time
  
    n <- sapply(1:ncomponents, function(j){sum(allocations==j)})
    temp = ds$weight - bGender*(as.numeric(ds$gender)-1) - bBy*(as.numeric(ds$by)-1)
    for(j in 1:nsubjects){
      startIndex = (j-1)*nrep + 1
      endIndex = j*nrep
      temp[startIndex:endIndex] = temp[startIndex:endIndex] - randComp[j,2]*time - randComp[j,1]
    }
    
    dgamma(errPrecision_max, shape = gammaShapeRate+nsubjects*nrep/2, 
           rate=gammaShapeRate+sum(temp^2)/2)
}))

######## Conditional posterior of beta #########
#p(beta|errPrecision, randPrecision, randmu, x)
mcmcfit_beta = fitModel_beta(niter, nthin, nburnin, jagsmodel = model_beta, 
                         nchains = numchains, ncomponents, randmu_max, 
                         randPrecision_max, errPrecision_max)$mcmcfit

chib = chib - log(mean(foreach(i=1:nrow(mcmcfit_beta[[1]]), .combine = c) %dopar%{
  
  ncomponents = attributes(mcmcfit_beta)$ncomponents
  nsubjects = attributes(mcmcfit_beta)$nsubjects
  time = attributes(mcmcfit_beta)$time
  nrep = attributes(mcmcfit_errPrecision)$nrep
  
  randComp = getRandomComp(mcmcfit_beta, i)
  allocations = getAllocations(mcmcfit_beta,i)
  
  temp = ds$weight
  for(j in 1:nsubjects){
    startIndex = (j-1)*nrep + 1
    endIndex = j*nrep
    temp[startIndex:endIndex] = temp[startIndex:endIndex] - randComp[j,2]*time - randComp[j,1]
  }
  
  tempreg = lm(temp~I(as.numeric(ds$gender)-1)+I(as.numeric(ds$by)-1)+0)
  designMatrix = model.matrix(tempreg)
  dmvnorm(c(betaGender_max, betaBy_max), mean = tempreg$coefficients, 
          sigma = solve(t(designMatrix)%*%designMatrix)/errPrecision_max)
}))

######## Conditional posterior of eta #########
#p(eta|errPrecision, randPrecision, randmu, beta, x)
mcmcfit_eta = fitModel_eta(niter, nthin, nburnin, jagsmodel = model_eta, 
                                             nchains = numchains, ncomponents, randmu_max, 
                                             randPrecision_max, errPrecision_max, betaGender_max,
                                            betaBy_max, betaTime_max)$mcmcfit

chib = chib - log(mean(foreach(i=1:nrow(mcmcfit_eta[[1]]), .combine = c) %dopar%{
  
  ncomponents = attributes(mcmcfit_eta)$ncomponents
  nsubjects = attributes(mcmcfit_eta)$nsubjects
  
  allocations = getAllocations(mcmcfit_eta,i)
  dirichParm=getDirichParam(ncomponents)
  nPerGroup <- sapply(1:ncomponents, function(j){sum(allocations==j)})
  
  exp(lgamma(sum(dirichParm)+nsubjects) - sum(lgamma(dirichParm+nPerGroup)) + 
        sum((dirichParm+nPerGroup-1)*log(eta_max)))
}))



