ppcCheck=foreach(k=1:nrow(mcmcfit[[1]]),.combine='c', .packages='MASS') %dopar%{
  
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
  betaAge = mcmcfit[[1]][k, "betaAge"]
  errSd = sqrt(1/mcmcfit[[1]][k, "errPrecision"])
  sigma = getSigma(mcmcfit, mcmcIterNum = k)
  
  #Only first data point
  ppRandomPart = matrix(nrow=0, ncol=nrep)
  
  totalPp = 0
  for(j in 1:ncomponents){
    numObs=round(eta[j]*10000)
    if(numObs>0){
      newObsRandPart = matrix(mvrnorm(n=numObs, randmu[[j]], randSigma[[j]]), ncol=2)
      
      for(m in 1:numObs){
        newObs = mvrnorm(n=1, newObsRandPart[m, 1] + newObsRandPart[m, 2]*time, diag(nrep)*errSd^2)
        newObs = newObs - mean(newObs)
        totalPp =totalPp + sum(newObs^2)
      }
    }
  }
  totalPp = totalPp/(10000*nrep)
  
  return(totalPp)
}

betaGender = mcmcfit[[1]][1, "betaGender"]
betaBy = mcmcfit[[1]][1, "betaBy"]
betaAge = mcmcfit[[1]][1, "betaAge"]
totalSample = c()
for(j in 1:nsubjects){
  fixedPartXBeta = betaBy * (as.numeric(dswide[j, "by"])-1) +
    betaGender * (as.numeric(dswide[j, "gender"])-1) + betaAge*dswide[j, "age"]
  
  sampleRandomPart = dswide_y[j,] - fixedPartXBeta
  sampleRandomPart = sampleRandomPart - mean(sampleRandomPart)
  totalSample[j] = sum(sampleRandomPart^2)
}
totalSample = sum(totalSample)/(nsubjects*nrep)

qplot(ppcCheck, geom=c("density"), xlab="Test statistic", ylab="PDF function estimated using KDE") + geom_vline(xintercept=totalSample, color="red")
plot(density(ppcCheck))
lines(density(ppcCheck[,1]), col='red')
