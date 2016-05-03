ppcCheck=foreach(k=1:nrow(mcmcfit[[1]]),.combine='c', .packages='mvtnorm') %dopar%{
  
  ncomponents = attributes(mcmcfit)$ncomponents
  nsubjects = attributes(mcmcfit)$nsubjects
  nrep = attributes(mcmcfit)$nrep
  time = attributes(mcmcfit)$time
  
  eta = getEta(mcmcfit, mcmcIterNum = k)
  randmu = getRandMu(mcmcfit, mcmcIterNum = k)
  randSigma = getRandSigma(mcmcfit, mcmcIterNum = k)
  allocations = getAllocations(mcmcfit, mcmcIterNum = k)
  randComp = getRandomComp(mcmcfit, mcmcIterNum = k)
  betaGender = mcmcfit[[1]][k, "betaGender"]
  betaBy = mcmcfit[[1]][k, "betaBy"]
  errSd = sqrt(1/mcmcfit[[1]][k, "errPrecision"])
  sigma = getSigma(mcmcfit, mcmcIterNum = k)
  
  #Only first data point
  tIndex = 10
  sampleFixedPart = ppFixedPart = matrix(nrow=nsubjects, ncol=nrep)
  total = 0
  for(j in 1:nsubjects){
    fixedPartXBeta = betaBy * (as.numeric(dswide[j, "by"])-1) + 
      betaGender * (as.numeric(dswide[j, "gender"])-1)
    
    randomPartXBeta = randmu[[allocations[j]]][1] + randmu[[allocations[j]]][2]*time
    #randomPartXBeta = randComp[j, 1] + randComp[j, 2]*time
    
    #op = rmvnorm(1, mean=fixedPartXBeta + randomPartXBeta, sigma = diag(nrep)*errSd^2)
    op = rmvnorm(1, mean=fixedPartXBeta + randomPartXBeta, sigma = sigma[[allocations[j]]])
    total = total + sum((op-mean(op))^2)
  }
  
  return(total/(nsubjects*nrep))
}

#calculating observed
total=0
for(j in 1:nsubjects){
  total = total + sum((dswide_y[j,]-mean(dswide_y[j,]))^2)
}
total = total/(nsubjects*nrep)

plot(density(ppcCheck))
abline(v = total)
HPDinterval(mcmc(ppcCheck))
