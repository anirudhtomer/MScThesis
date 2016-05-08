ppcCheck=foreach(k=1:nrow(mcmcfit[[1]]),.combine='c', .packages='MASS') %dopar%{
  
  ncomponents = attributes(mcmcfit)$ncomponents
  nsubjects = attributes(mcmcfit)$nsubjects
  INTERCEPT_SCALE = attributes(mcmcfit)$INTERCEPT_SCALE

  eta = getEta(mcmcfit, mcmcIterNum = k)
  randmu = getRandMu(mcmcfit, mcmcIterNum = k)
  randSigma = getRandSigma(mcmcfit, mcmcIterNum = k)
  allocations = getAllocations(mcmcfit, mcmcIterNum = k)
  randComp = getRandomComp(mcmcfit, mcmcIterNum = k)
  errSd = sqrt(1/mcmcfit[[1]][k, "errPrecision"])
  
  betaAge = mcmcfit[[1]][k,"betaAge"]
  betaDonate = mcmcfit[[1]][k,"betaDonate"]
  betaSeason = mcmcfit[[1]][k,"betaSeason"]
  betaTSPD = mcmcfit[[1]][k,"betaTSPD"]
  betaDonateLast2TSPD = mcmcfit[[1]][k,"betaDonateLast2TSPD"]
  betaDonateLast2Donate = mcmcfit[[1]][k,"betaDonateLast2Donate"]
  betaDonateLast2Square = mcmcfit[[1]][k,"betaDonateLast2Square"]
  
  nrep = 15
  totalPp = 0
  for(j in 1:ncomponents){
    numObs=round(eta[j]*10000)
    if(numObs>0){
      newObsRandPart = matrix(mvrnorm(n=numObs, randmu[[j]], randSigma[[j]]), ncol=2)
      
      for(m in 1:numObs){
        donationLast2years = round(rnorm(nrep, mean = sample(1:10, size = 1), sd = 1))
        donationLast2years[donationLast2years<0] = 0
        donationLast2years = donationLast2years/100
        newObs = mvrnorm(n=1, newObsRandPart[m, 1]*INTERCEPT_SCALE + newObsRandPart[m, 2]*donationLast2years, diag(nrep)*errSd^2)
        newObs = newObs - mean(newObs)
        totalPp =totalPp + sum(newObs^2)
      }
    }
  }
  totalPp = totalPp/(10000*nrep)
  
  return(totalPp)
}



temp = c()
for(i in 1:nrow(mcmcfit[[1]])){
randComp = getRandomComp(mcmcfit, mcmcIterNum = i)

betaAge = mcmcfit[[1]][i,"betaAge"]
betaDonate = mcmcfit[[1]][i,"betaDonate"]
betaSeason = mcmcfit[[1]][i,"betaSeason"]
betaTSPD = mcmcfit[[1]][i,"betaTSPD"]
betaDonateLast2TSPD = mcmcfit[[1]][i,"betaDonateLast2TSPD"]
betaDonateLast2Donate = mcmcfit[[1]][i,"betaDonateLast2Donate"]
betaDonateLast2Square = mcmcfit[[1]][i,"betaDonateLast2Square"]

betaAge = -0.08576
betaDonate = 0.1546
betaSeason = -0.08014
betaTSPD =-0.03641
betaDonateLast2TSPD = 0.01964
betaDonateLast2Donate = 4.6887
betaDonateLast2Square = -52.7956

totalSample = c()
for(j in 1:nsubjects){
  startIndex = cumsumHb[j]+1
  endIndex = cumsumHb[j]+numHb[j]
  
  fixedPartMean = betaAge*ds$Age[startIndex:endIndex] + betaSeason*(as.numeric(ds$Season[startIndex:endIndex])-1) +
    betaDonate*(as.numeric(ds$Donate[startIndex:endIndex])) + betaTSPD*ds$TSPD[startIndex:endIndex] +
    betaDonateLast2TSPD*ds$donateLast2TSPD[startIndex:endIndex] +
    betaDonateLast2Donate*ds$donateLast2Donate[startIndex:endIndex] +
    betaDonateLast2Square*ds$donateLast2Square[startIndex:endIndex]
  
  sampleRandomPart = ds$Hb[startIndex:endIndex] - fixedPartMean
  sampleRandomPart = sampleRandomPart - mean(sampleRandomPart)
  totalSample[j] = sum(sampleRandomPart^2)
}
totalSample = sum(totalSample)/sum(numHb)
temp[i] = totalSample
}

qplot(ppcCheck_5comp, geom=c("density"), xlab="Test statistic", ylab="PDF function estimated using KDE") + geom_vline(xintercept=totalSample, color="red")

plot(density(ppcCheck))
lines(density(temp))
