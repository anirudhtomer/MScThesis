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
randPrecision_max = getRandPrecision(mcmcfit, mcmcIterNum = maxValueIter)

fixedEffects_max = c("betaGender"= mcmcfit[[1]][maxValueIter, "betaGender"], 
                     "betaBy"=mcmcfit[[1]][maxValueIter, "betaBy"])
errPrecision_max = mcmcfit[[1]][maxValueIter, "errPrecision"]
eta_max = getEta(mcmcfit, mcmcIterNum = maxValueIter)

#Step 2: log prior p(theta*)
log_prior=lgamma(sum(dirichParm)) - sum(lgamma(dirichParm)) + sum((dirichParm-1)*log(eta_max)) + 
  sum(dnorm(c(randIntercept_max, randSlope_max, fixedEffects_max), 
            mean = betaMu, sd = sqrt(1/betaTau), log = T)) + 
  dgamma(errPrecision_max, shape = gammaShapeRate, rate = gammaShapeRate, log = T) + 
  sum(sapply(randPrecision_max, function(x){log(dwish(x, v = wishartPriorDf, S = wishartPriorScale))}))

chib = obslogL[maxValueIter] + log_prior

###### Conditional posterior of randprecision #######
#p(randPrecision|y) = E(p(randPrecision|z,y,mu)) over p(z,mu|y)
mcmcfit_randPrecision = mcmcfit
chib = chib - log(mean(foreach(i=1:nrow(mcmcfit[[1]]), .combine = c) %dopar%{
  
  allocations = getAllocations(mcmcfit_randPrecision, mcmcIterNum = i)
  randComp = getRandomComp(mcmcfit_randPrecision, mcmcIterNum = i)
  ncomponents = attributes(mcmcfit_randPrecision)$ncomponents
  randmu = getRandMu(mcmcfit_randPrecision, mcmcIterNum = i)
  
  freqPerGroup = sapply(1:ncomponents, function(x){sum(allocations==x)})
  randCompPerGroup = lapply(1:ncomponents, function(x){randComp[allocations==x,]})
  
  wishartPostDf = freqPerGroup + wishartPriorDf
  wishartPosteriorScale = rep(list(wishartPriorScale), ncomponents)
  
  for(j in 1:ncomponents){
    if(freqPerGroup[j]>0){
      x_minus_mu = apply(randCompPerGroup[[j]], MARGIN=1, FUN=function(x){x-randmu[[j]]})
      wishartPosteriorScale[[j]] = solve(solve(wishartPriorScale) + x_minus_mu%*%t(x_minus_mu))
    }
  }
  
  prob = 1
  for(j in 1:ncomponents){
    prob = prob * dwish(randPrecision_max[[j]], v = wishartPostDf[j], S = wishartPosteriorScale[[j]])
  }
  
  return(prob)
}))

######## Conditional posterior of randmu #########
#p(randMu|randPrecision_max,y) = p(randmu|z,y, randPrecision_max) p(z|y, randPrecision_max)
mcmcfit_randmu = fitModel_randmu(niter, nthin, nburnin, jagsmodel = model_randmu, nchains = numchains, ncomponents,
                                 randPrecision_max)$mcmcfit

chib = chib - log(mean(foreach(i=1:nrow(mcmcfit_randmu[[1]]), .combine = c) %dopar%{

  allocations = getAllocations(mcmcfit_randmu, mcmcIterNum = i)
  randComp = getRandomComp(mcmcfit_randmu, mcmcIterNum = i)
  ncomponents = attributes(mcmcfit_randmu)$ncomponents
  randmu = getRandMu(mcmcfit_randmu, mcmcIterNum = i)
  
  freqPerGroup = sapply(1:ncomponents, function(x){sum(allocations==x)})
  randCompPerGroup = lapply(1:ncomponents, function(x){randComp[allocations==x,]})
  
  prob = 1
  for(j in 1:ncomponents){
    sigmaPost = solve(diag(2)*betaTau + freqPerGroup[j]*randPrecision_max[[j]])
    
    meanPost = c(betaMu, betaMu)
    if(freqPerGroup[j]>0){
      meanPost = c(sigmaPost%*%(
          (diag(2)*betaTau)%*%meanPost + freqPerGroup[j]*randPrecision_max[[j]]%*%apply(randCompPerGroup[[j]], MARGIN = 2, FUN = mean)
          ))
    }
    prob = prob * dmvnorm(x = c(randIntercept_max, randSlope_max), mean = meanPost, sigma = sigmaPost,log = F)
  }
  
  return(prob)
}))

######## Conditional posterior of errPrecision #########
#p(errPrecision|randPrecision_max, randmu_max, y) = p(errPrecision|z,y,beta,randPrecision_max, randmu_max) p(z,beta|y,randPrecision_max, randmu_max)
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
#p(beta|errPrecision_max, randPrecision_max, randmu_max, y) = p(beta|z,y,errPrecision_max,randPrecision_max, randmu_max) p(z|y, beta_max, randPrecision_max, randmu_max)
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



