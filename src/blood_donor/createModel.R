model = function(){
  for(i in 1:nsubjects){
    for(j in 1:numHb[i]){
      dsHb[cumsumHb[i]+j]~dnorm(mu[cumsumHb[i]+j], errPrecision)
    
      mu[cumsumHb[i]+j]<- betaDonate*dsDonate[cumsumHb[i]+j] + betaAge*dsAge[cumsumHb[i]+j] + 
        betaSeason*dsSeason[cumsumHb[i]+j] + betaTSPD*dsTSPD[cumsumHb[i]+j] +
        randomComp[i,1]*0.1 + randomComp[i,2]*(dsDonationLast2Years[cumsumHb[i]+j]) + 
        betaAgeSeason*dsAgeSeason[cumsumHb[i]+j] + 
        betaTSPDDonate*dsTSPDDonate[cumsumHb[i]+j] +
        betaTSPDSeason*dsTSPDSeason[cumsumHb[i]+j] +
        
        betaDonateLast2TSPD*dsDonateLast2TSPD[cumsumHb[i]+j] +
        betaDonateLast2Donate*dsDonateLast2Donate[cumsumHb[i]+j] +
        betaDonateLast2Season*dsDonateLast2Season[cumsumHb[i]+j] +
        betaDonateLast2Square*dsDonateLast2Square[cumsumHb[i]+j]
    }
    
    randomComp[i,1:2]~dmnorm(randmu[S[i], 1:2], randPrecision[1:2,1:2,S[i]])
    S[i]~dcat(Eta[])
  }
  
  for(k in 1:ncomponents){
    randmu_unordered[k,1]~dnorm(betaMu, betaTau)
    randmu_unordered[k,2]~dnorm(betaMu, betaTau)
    
    randPrecision[1:2,1:2,k]~dwish(wishartPriorScale, wishartPriorDf)
    randSigma[1:2,1:2,k]<-inverse(randPrecision[1:2,1:2,k])
    
    precisionIntercept[k]<-randPrecision[1,1,k]
    precisionSlope[k]<-randPrecision[2,2,k]
    rho[k]<-randSigma[1,2,k]* sqrt(precisionIntercept[k] * precisionSlope[k])
  }
  
  orderRandMu<-order(randmu_unordered[1:ncomponents,1])
  for(h in 1:ncomponents){
    randmu[h,1:2]<-randmu_unordered[orderRandMu[h],1:2]
  }
  
  errPrecision~dgamma(gammaShapeRate,gammaShapeRate)
  
  betaDonate~dnorm(betaMu,betaTau)
  betaAge~dnorm(betaMu,betaTau)
  betaSeason~dnorm(betaMu,betaTau)
  betaTSPD~dnorm(betaMu,betaTau)
  
  betaAgeSeason~dnorm(betaMu,betaTau)
  betaTSPDDonate~dnorm(betaMu,betaTau)
  betaTSPDSeason~dnorm(betaMu,betaTau)
  betaDonateLast2TSPD~dnorm(betaMu,betaTau)
  betaDonateLast2Donate~dnorm(betaMu,betaTau)
  betaDonateLast2Season~dnorm(betaMu,betaTau)
  betaDonateLast2Square~dnorm(betaMu,betaTau)
  
  Eta~ddirch(dirichParm[])
}

singleModel = function(){
  for(i in 1:nsubjects){
    for(j in 1:numHb[i]){
      dsHb[cumsumHb[i]+j]~dnorm(mu[cumsumHb[i]+j], errPrecision)
      
      mu[cumsumHb[i]+j]<- betaDonate*dsDonate[cumsumHb[i]+j] + betaAge*dsAge[cumsumHb[i]+j] + 
        betaSeason*dsSeason[cumsumHb[i]+j] + betaTSPD*dsTSPD[cumsumHb[i]+j] +
        randomComp[i,1]*0.1 + randomComp[i,2]*(dsDonationLast2Years[cumsumHb[i]+j]) + 
        
        betaAgeSeason*dsAgeSeason[cumsumHb[i]+j] + 
        betaTSPDDonate*dsTSPDDonate[cumsumHb[i]+j] +
        betaTSPDSeason*dsTSPDSeason[cumsumHb[i]+j] +
        
        betaDonateLast2TSPD*dsDonateLast2TSPD[cumsumHb[i]+j] +
        betaDonateLast2Donate*dsDonateLast2Donate[cumsumHb[i]+j] +
        betaDonateLast2Season*dsDonateLast2Season[cumsumHb[i]+j] +
        betaDonateLast2Square*dsDonateLast2Square[cumsumHb[i]+j]
    }
    
    randomComp[i,1:2]~dmnorm(randmu[S[i], 1:2], randPrecision[1:2,1:2,S[i]])
    S[i]~dcat(Eta[])
  }
  
  randmu[1,1]~dnorm(betaMu, betaTau)
  randmu[1,2]~dnorm(betaMu, betaTau)
  
  randPrecision[1:2,1:2,1]~dwish(wishartPriorScale, wishartPriorDf)
  randSigma[1:2,1:2,1]<-inverse(randPrecision[1:2,1:2,1])
  
  precisionIntercept[1]<-randPrecision[1,1,1]
  precisionSlope[1]<-randPrecision[2,2,1]
  rho[1]<-randSigma[1,2,1]* sqrt(precisionIntercept[1] * precisionSlope[1])
  
  errPrecision~dgamma(gammaShapeRate,gammaShapeRate)
  
  betaDonate~dnorm(betaMu,betaTau)
  betaAge~dnorm(betaMu,betaTau)
  betaSeason~dnorm(betaMu,betaTau)
  betaTSPD~dnorm(betaMu,betaTau)
  
  betaAgeSeason~dnorm(betaMu,betaTau)
  betaTSPDDonate~dnorm(betaMu,betaTau)
  betaTSPDSeason~dnorm(betaMu,betaTau)
  betaDonateLast2TSPD~dnorm(betaMu,betaTau)
  betaDonateLast2Donate~dnorm(betaMu,betaTau)
  betaDonateLast2Season~dnorm(betaMu,betaTau)
  betaDonateLast2Square~dnorm(betaMu,betaTau)
  
  Eta[1]~dbeta(100,1)
}

