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
  
  bloodDonorId = unique(dsall$Id)
  
  nrep = 15
  totalPp = 0
  totalEntries = 0
  for(j in 1:ncomponents){
    numObs=round(eta[j]*10000)
    if(numObs>0){
      newObsRandPart = matrix(mvrnorm(n=numObs, randmu[[j]], randSigma[[j]]), ncol=2)
      
      for(m in 1:numObs){
        tempId = bloodDonorId[sample((1:length(bloodDonorId)), size = 1)]
        tempdonationsLast2years = dsall[dsall$Id == tempId, "donationLast2Years"]
        nrep = length(tempdonationsLast2years)
        totalEntries = totalEntries + nrep
        
        newObs = mvrnorm(n=1, newObsRandPart[m, 1]*INTERCEPT_SCALE + newObsRandPart[m, 2]*tempdonationsLast2years, diag(nrep)*errSd^2)
        newObs = newObs - mean(newObs)
        totalPp = totalPp + sum(newObs^2)
      }
    }
  }
  totalPp = totalPp/totalEntries
  
  return(totalPp)
}

temp = c()
for(i in 1:nrow(mcmcfit[[1]])){
randComp = getRandomComp(mcmcfit, mcmcIterNum = i)

betaAge = mcmcfit[[1]][i,"betaAge"]
betaDonate = mcmcfit[[1]][i,"betaDonate"]
betaSeason = mcmcfit[[1]][i,"betaSeason"]
betaTSPD = mcmcfit[[1]][i,"betaTSPD"]
# betaTSPDSquare = mcmcfit[[1]][i,"betaTSPDSquare"]
# betaTSPDSeason = mcmcfit[[1]][i,"betaTSPDSeason"]
betaDonateLast2TSPD = mcmcfit[[1]][i,"betaDonateLast2TSPD"]
betaDonateLast2Donate = mcmcfit[[1]][i,"betaDonateLast2Donate"]
betaDonateLast2Square = mcmcfit[[1]][i,"betaDonateLast2Square"]

totalSample = c()
for(j in 1:nsubjects){
  startIndex = cumsumHb[j]+1
  endIndex = cumsumHb[j]+numHb[j]
  
  fixedPartMean = betaAge*ds$Age[startIndex:endIndex] + betaSeason*(as.numeric(ds$Season[startIndex:endIndex])-1) +
    betaDonate*(as.numeric(ds$Donate[startIndex:endIndex])) + betaTSPD*ds$TSPD[startIndex:endIndex] +
    betaDonateLast2TSPD*ds$donateLast2TSPD[startIndex:endIndex] +
    betaDonateLast2Donate*ds$donateLast2Donate[startIndex:endIndex] +
    betaDonateLast2Square*ds$donateLast2Square[startIndex:endIndex]
    # betaTSPDSquare*ds$TSPD[startIndex:endIndex]*ds$TSPD[startIndex:endIndex] +
    # betaTSPDSeason*ds$TSPDSeason[startIndex:endIndex]
  
  sampleRandomPart = ds$Hb[startIndex:endIndex] - fixedPartMean
  sampleRandomPart = sampleRandomPart - mean(sampleRandomPart)
  totalSample[j] = sum(sampleRandomPart^2)
}
totalSample = sum(totalSample)/sum(numHb)
temp[i] = totalSample
}

resultDf = data.frame("Test.statistic"=c(ppcCheck, temp), "Type"=c(rep("Posterior Predictive", length(ppcCheck)), rep("Sample", length(temp))))

qplot(Test.statistic, geom=c("density"), data=resultDf, xlab="Test statistic", ylab="PDF function estimated using KDE", color=Type) + 
  scale_x_continuous(breaks = seq(0, 3, by = 0.25))
