#Based on Chib's approximation
source("DIC_functions.R")

copyAttributesAndTrimLength=function(destination, source){
  attributes(destination)$ncomponents = attributes(source)$ncomponents
  attributes(destination)$nsubjects = attributes(source)$nsubjects
  attributes(destination)$nrep = attributes(source)$nrep
  attributes(destination)$time = attributes(source)$time
  attributes(destination)$wishartPriorScale = attributes(source)$wishartPriorScale
  if(attributes(destination)$ncomponents==1){
    colnames(destination[[1]])[match("Eta",colnames(destination[[1]]))] = "Eta[1]"
  }
  
  destination[[1]] = destination[[1]][((nburnin/nthin+1):(niter/nthin)),]
  attributes(destination[[1]])$mcpar = c(1, nrow(destination[[1]]), 1)
  return(destination)
}

#Step 1: logL
#calculated via the DIC function calculateObsTotalDTheta
#since the function I am calling gives -2*logL, I multiply it with -0.5
obslogL=calculateObsTotalDTheta(mcmcfit)*(-0.5)
maxValueIter = which.max(obslogL)

randmu_max = getRandMu(mcmcfit, mcmcIterNum = maxValueIter)
randIntercept_max = sapply(randmu_max, function(x){x[1]})
randSlope_max = sapply(randmu_max, function(x){x[2]})
randPrecision_max = getRandPrecision(mcmcfit, mcmcIterNum = maxValueIter)
randSigma_max = getRandSigma(mcmcfit, mcmcIterNum = maxValueIter)

betaGender_max = mcmcfit[[1]][maxValueIter, "betaGender"]
betaBy_max = mcmcfit[[1]][maxValueIter, "betaBy"]
betaAge_max = mcmcfit[[1]][maxValueIter, "betaAge"]

errPrecision_max = mcmcfit[[1]][maxValueIter, "errPrecision"]
eta_max = getEta(mcmcfit, mcmcIterNum = maxValueIter)

dirichParm = rep(1, attributes(mcmcfit)$ncomponents)
#Step 2: log prior p(theta*)
log_prior=lgamma(sum(dirichParm)) - sum(lgamma(dirichParm)) + sum((dirichParm-1)*log(eta_max)) + 
  sum(dnorm(c(randIntercept_max, randSlope_max, betaGender_max, betaBy_max, betaAge_max), 
            mean = betaMu, sd = sqrt(1/betaTau), log = T)) + 
  dgamma(errPrecision_max, shape = gammaShapeRate, rate = gammaShapeRate, log = T) + 
  sum(sapply(randPrecision_max, function(x){log(dwish(x, v = wishartPriorDf, S = attributes(mcmcfit)$wishartPriorScale))}))

chib = obslogL[maxValueIter] + log_prior

###### Conditional posterior of randprecision #######
#p(randPrecision|y) = E(p(randPrecision|z,y,mu)) over p(z,mu|y)
mcmcfit_randPrecision = mcmcfit
chib = chib - log(mean(foreach(i=1:nrow(mcmcfit[[1]]), .combine = c, .packages="MCMCpack") %dopar%{
  
  allocations = getAllocations(mcmcfit_randPrecision, mcmcIterNum = i)
  randComp = getRandomComp(mcmcfit_randPrecision, mcmcIterNum = i)
  ncomponents = attributes(mcmcfit_randPrecision)$ncomponents
  randmu = getRandMu(mcmcfit_randPrecision, mcmcIterNum = i)
  
  freqPerGroup = sapply(1:ncomponents, function(x){sum(allocations==x)})
  randCompPerGroup = lapply(1:ncomponents, function(x){matrix(randComp[allocations==x,], ncol=2)})
  
  wishartPostDf = freqPerGroup + wishartPriorDf
  wishartPosteriorScale = rep(list(attributes(mcmcfit_randPrecision)$wishartPriorScale), ncomponents)
  
  for(j in 1:ncomponents){
    if(freqPerGroup[j]>0){
      x_minus_mu = apply(randCompPerGroup[[j]], MARGIN=1, FUN=function(x){x-randmu[[j]]})
      wishartPosteriorScale[[j]] = solve(solve(attributes(mcmcfit_randPrecision)$wishartPriorScale) + x_minus_mu%*%t(x_minus_mu))
    }
  }
  
  prob = 1
  for(j in 1:ncomponents){
    prob = prob * dWISHART(randPrecision_max[[j]], df = wishartPostDf[j], S = wishartPosteriorScale[[j]])
  }
  
  return(prob)
}))

######## Conditional posterior of randmu #########
#p(randMu|randPrecision_max,y) = p(randmu|z,y, randPrecision_max) p(z|y, randPrecision_max)
if(attributes(mcmcfit)$ncomponents==1){
  fit_randmu = fitModel_randmu(niter, nthin, 0, jagsmodel = singleModel_randmu, nchains = numchains, attributes(mcmcfit)$ncomponents,
                               randPrecision_max)
  
}else{
  fit_randmu = fitModel_randmu(niter, nthin, 0, jagsmodel = model_randmu, nchains = numchains, attributes(mcmcfit)$ncomponents,
                               randPrecision_max)
}
mcmcfit_randmu = as.mcmc(fit_randmu)
mcmcfit_randmu = copyAttributesAndTrimLength(mcmcfit_randmu, mcmcfit)

chib = chib - log(mean(foreach(i=1:nrow(mcmcfit_randmu[[1]]), .combine = c) %dopar%{

  allocations = getAllocations(mcmcfit_randmu, mcmcIterNum = i)
  randComp = getRandomComp(mcmcfit_randmu, mcmcIterNum = i)
  ncomponents = attributes(mcmcfit_randmu)$ncomponents
  randmu = getRandMu(mcmcfit_randmu, mcmcIterNum = i)
  
  freqPerGroup = sapply(1:ncomponents, function(x){sum(allocations==x)})
  randCompPerGroup = lapply(1:ncomponents, function(x){matrix(randComp[allocations==x,], ncol=2)})
  
  prob = 1
  for(j in 1:ncomponents){
    sigmaPost = solve(diag(2)*betaTau + freqPerGroup[j]*randPrecision_max[[j]])
    
    meanPost = c(betaMu, betaMu)
    if(freqPerGroup[j]>0){
      meanPost = c(sigmaPost%*%(
          (diag(2)*betaTau)%*%meanPost + freqPerGroup[j]*randPrecision_max[[j]]%*%apply(randCompPerGroup[[j]], MARGIN = 2, FUN = mean)
          ))
    }
    
    prob = prob * dmvnorm(x = c(randIntercept_max[j], randSlope_max[j]), mean = meanPost, sigma = sigmaPost,log = F)
  }
  
  return(prob)
}))

######## Conditional posterior of errPrecision #########
#p(errPrecision|randPrecision_max, randmu_max, y) = p(errPrecision|z,y,beta,randPrecision_max, randmu_max) p(z,beta|y,randPrecision_max, randmu_max)
if(attributes(mcmcfit)$ncomponents==1){
  fit_errPrecision = fitModel_errPrecision(niter, nthin, 0, jagsmodel = singleModel_errPrecision,
                        nchains = numchains, attributes(mcmcfit)$ncomponents, randmu_max, randPrecision_max)
}else{
  fit_errPrecision = fitModel_errPrecision(niter, nthin, 0, jagsmodel = model_errPrecision,
                        nchains = numchains, attributes(mcmcfit)$ncomponents, randmu_max, randPrecision_max)
  
}
mcmcfit_errPrecision = as.mcmc(fit_errPrecision)
mcmcfit_errPrecision = copyAttributesAndTrimLength(mcmcfit_errPrecision, mcmcfit)

chib = chib - log(mean(foreach(i=1:nrow(mcmcfit_errPrecision[[1]]), .combine = c) %dopar%{

    allocations = getAllocations(mcmcfit_randmu, mcmcIterNum = i)
    randComp = getRandomComp(mcmcfit_randmu, mcmcIterNum = i)
    bGender = mcmcfit_errPrecision[[1]][i, "betaGender"]
    bBy = mcmcfit_errPrecision[[1]][i, "betaBy"]
    bAge = mcmcfit_errPrecision[[1]][i, "betaAge"]
    
    ncomponents = attributes(mcmcfit_errPrecision)$ncomponents
    nsubjects = attributes(mcmcfit_errPrecision)$nsubjects
    nrep = attributes(mcmcfit_errPrecision)$nrep
    time = attributes(mcmcfit_errPrecision)$time
  
    temp = ds$weight - bGender*(as.numeric(ds$gender)-1) - bBy*(as.numeric(ds$by)-1) - bAge*ds$age
    for(j in 1:nsubjects){
      startIndex = (j-1)*nrep + 1
      endIndex = j*nrep
      temp[startIndex:endIndex] = temp[startIndex:endIndex] - randComp[j,2]*time - randComp[j,1]
    }
    
    dgamma(errPrecision_max, shape = gammaShapeRate+nsubjects*nrep/2, 
           rate=gammaShapeRate+sum(temp^2)/2)
}))

######## Conditional posterior of beta #########
#p(beta|errPrecision_max, randPrecision_max, randmu_max, y) = p(beta|z,y,errPrecision_max,randPrecision_max, randmu_max) p(z|y, beta_max, randPrecision_max, randmu_max)
if(attributes(mcmcfit)$ncomponents==1){
  fit_beta =  fitModel_beta(niter, nthin, 0, jagsmodel = singleModel_beta, 
                                    nchains = numchains, attributes(mcmcfit)$ncomponents, randmu_max, 
                                    randPrecision_max, errPrecision_max)
}else{
  fit_beta =  fitModel_beta(niter, nthin, 0, jagsmodel = model_beta, 
                                    nchains = numchains, attributes(mcmcfit)$ncomponents, randmu_max, 
                                    randPrecision_max, errPrecision_max)
  
}
mcmcfit_beta = as.mcmc(fit_beta)
mcmcfit_beta = copyAttributesAndTrimLength(mcmcfit_beta, mcmcfit)

chib = chib - log(mean(foreach(i=1:nrow(mcmcfit_beta[[1]]), .combine = c) %dopar%{
  
  ncomponents = attributes(mcmcfit_beta)$ncomponents
  nsubjects = attributes(mcmcfit_beta)$nsubjects
  time = attributes(mcmcfit_beta)$time
  nrep = attributes(mcmcfit_beta)$nrep
  
  randComp = getRandomComp(mcmcfit_beta, mcmcIterNum=i)
  allocations = getAllocations(mcmcfit_beta,mcmcIterNum=i)
  
  temp = ds$weight
  for(j in 1:nsubjects){
    startIndex = (j-1)*nrep + 1
    endIndex = j*nrep
    temp[startIndex:endIndex] = temp[startIndex:endIndex] - randComp[j,2]*time - randComp[j,1]
  }
  
  tempreg = lm(temp~I(as.numeric(ds$gender)-1)+I(as.numeric(ds$by)-1)+ds$age+0)
  designMatrix = model.matrix(tempreg)
  dmvnorm(c(betaGender_max, betaBy_max, betaAge_max), mean = tempreg$coefficients, 
          sigma = solve(t(designMatrix)%*%designMatrix)/errPrecision_max)
}))

######## Conditional posterior of eta #########
#p(eta|errPrecision, randPrecision, randmu, beta, y)
if(attributes(mcmcfit)$ncomponents==1){
  fit_eta =  fitModel_eta(niter, nthin, 0, jagsmodel = singleModel_eta, 
                          nchains = numchains, attributes(mcmcfit)$ncomponents, randmu_max, 
                          randPrecision_max, errPrecision_max, betaGender_max,
                          betaBy_max, betaAge_max)
}else{
  fit_eta =  fitModel_eta(niter, nthin, 0, jagsmodel = model_eta, 
                           nchains = numchains, attributes(mcmcfit)$ncomponents, randmu_max, 
                           randPrecision_max, errPrecision_max, betaGender_max,
                           betaBy_max, betaAge_max)
  
}
mcmcfit_eta = as.mcmc(fit_eta)
mcmcfit_eta = copyAttributesAndTrimLength(mcmcfit_eta, mcmcfit)

chib = chib - log(mean(foreach(i=1:nrow(mcmcfit_eta[[1]]), .combine = c) %dopar%{
  
  ncomponents = attributes(mcmcfit_eta)$ncomponents
  nsubjects = attributes(mcmcfit_eta)$nsubjects
  
  allocations = getAllocations(mcmcfit_beta,mcmcIterNum=i)
  nPerGroup <- sapply(1:ncomponents, function(j){sum(allocations==j)})
  
  exp(lgamma(sum(dirichParm)+nsubjects) - sum(lgamma(dirichParm+nPerGroup)) + 
        sum((dirichParm+nPerGroup-1)*log(eta_max)))
}))



