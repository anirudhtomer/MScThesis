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
  
  #Only first data point
  tIndex = 10
  sampleFixedPart = ppFixedPart = matrix(nrow=nsubjects, ncol=nrep)
  total = numeric()
  for(j in 1:nsubjects){
    fixedPartXBeta = betaBy * (as.numeric(dswide[i, "by"])-1) + 
      betaGender * (as.numeric(dswide[i, "gender"])-1) + betaAge*dswide[i, "age"]
    
    ppFixedPart[j, ] = fixedPartXBeta + rnorm(n = nrep, 0, sd = errSd)
   sampleFixedPart[j,] = dswide_y[j,] - (randComp[j, 1] + randComp[j, 2]*time)
   total = c(total, ppFixedPart[j,]-sampleFixedPart[j,])
  }
  
  # #only first data point i.e. t=1
  # ppRandomPart = numeric()
  # for(j in 1:ncomponents){
  #   tempNsubjects = round(eta[j]*10000)
  #   if(tempNsubjects > 0){
  #     temp = mvrnorm(n=tempNsubjects, 
  #                    mu = randmu[[j]], Sigma = randSigma[[j]])
  #     if(tempNsubjects==1){
  #       temp = matrix(temp, ncol=2)
  #     }
  #     randomPartZb = apply(temp, MARGIN = 1, function(rowData){rowData[1]+time[tIndex]*rowData[2]})
  #     ppRandomPart = c(ppRandomPart, randomPartZb+rnorm(length(randomPartZb), mean = 0, sd = errSd))
  #   }
  # }
  
  return(mean(total))
}
plot(density(ppcCheck))
HPDinterval(mcmc(ppcCheck))
