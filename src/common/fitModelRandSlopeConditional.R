initS=function(nfolks, ncomponents){
  retval = numeric(nfolks)
  randomorder=sample(1:nfolks, nfolks)
  for(i in 1:ncomponents){
    retval[randomorder[((i-1)*(nfolks/ncomponents)+1):(i*(nfolks/ncomponents))]] = i;
  }
  retval[retval==0]=1
  return (retval)
}

betaMu = 0
betaTau = 0.0001
gammaShapeRate = 0.0001
wishartPriorScale = diag(2)*10
wishartPriorDf = 3

runJagsModel = function(datanodes, initialValues, stochasticNodes,jagsmodel, ncomp, niter, nthin, nburnin, nchains=1){
  if(ncomp==1){
    #The following line removes the element Eta from initial values
    initialValues$Eta = NULL
    datanodes$dirichParm = NULL
    datanodes$S = initialValues$S
    initialValues$S=NULL
    datanodes$Eta=c(1)
    datanodes$ncomponents = NULL
  }
  
  unload.module("glm")
  chainInits = list(initialValues)
  nchainsleft = nchains-1
  
  while(nchainsleft > 0){
    chainInits[[nchains - nchainsleft + 1]] = initialValues
    nchainsleft = nchainsleft - 1
  }
  
  fit = jags(data=datanodes, inits=chainInits, stochasticNodes,
             n.chains=nchains, n.iter=niter, n.thin=nthin, n.burnin=nburnin,
             model.file=jagsmodel, jags.module=NULL)
}

fitModel = function(niter=10000, nthin=50, nburnin=200, nchains=1, jagsmodel, ncomp=3){
  
  dirichParm = rep(0.1, ncomp)
  randComp=extractRandomComp(viaReg = T)
  quantiles = seq(1/(ncomp+1),ncomp/(ncomp+1),length.out = ncomp)
  randIntQuantiles = quantile(randComp[,1], probs = quantiles)
  randSlopeQuantiles = quantile(randComp[,2], probs = quantiles)
  reg=lm(weight~gender+by+age+time,data=ds)
  
  datanodes = list("dsweight"=ds$weight,"dsby"=as.numeric(ds$by)-1, "dsage"=ds$age,
                   "dsgender"=as.numeric(ds$gender)-1,"dstime"=ds$time,
                   "nsubjects"=nsubjects,"nrep"=nrep, "gammaShapeRate"=gammaShapeRate, "betaMu"=betaMu,
                   "betaTau"=betaTau, "betaMuMult"=rep(betaMu,2),"ncomponents"=ncomp, "betaTauMult"=diag(2),
                  "dirichParm"=dirichParm, "wishartPriorScale"=wishartPriorScale, "wishartPriorDf"=wishartPriorDf)
  
  initialValues = list("betaGender"=c(reg$coefficients[2]),"betaBy"=c(reg$coefficients[3]), "betaAge"=c(reg$coefficients[4]),
                       "errPrecision"=c(1), 
                       #"precisionIntercept"=rep(1, ncomp), "precisionSlope"=rep(1, ncomp),"rho"=rep(0, ncomp),
                       "Eta"=rep(1/ncomp, ncomp),
                       "S"=initS(nsubjects, ncomp),
                       "randmu_unordered"=cbind(randIntQuantiles, randSlopeQuantiles))
  
  stochasticNodes = c("betaGender","betaBy", "betaAge", "errPrecision", "randSigma","randPrecision",
                      "randmu", "randomComp", "Eta", "S", "precisionIntercept","precisionSlope", "rho")
  
  
  return(runJagsModel(datanodes, initialValues, stochasticNodes, jagsmodel, ncomp, niter, nthin, nburnin, nchains))
}

fitModel_randmu = function(niter=10000, nthin=50, nburnin=200, nchains=1, jagsmodel, ncomp=3, randPrecision_max){
  
  dirichParm = rep(1, ncomp)
  randComp=extractRandomComp(viaReg = T)
  quantiles = seq(1/(ncomp+1),ncomp/(ncomp+1),length.out = ncomp)
  randIntQuantiles = quantile(randComp[,1], probs = quantiles)
  randSlopeQuantiles = quantile(randComp[,2], probs = quantiles)
  reg=lm(weight~gender+by+age+time,data=ds)
  
  randPrecision = array(data=NA, dim = c(2,2,ncomp))
  for(i in 1:ncomp){
    randPrecision[,,i] = randPrecision_max[[i]]
  }
  
  datanodes = list("dsweight"=ds$weight,"dsby"=as.numeric(ds$by)-1, "dsage"=ds$age,
                   "dsgender"=as.numeric(ds$gender)-1,"dstime"=ds$time,
                   "nsubjects"=nsubjects,"nrep"=nrep, "gammaShapeRate"=gammaShapeRate, "betaMu"=betaMu,
                   "betaTau"=betaTau, "betaMuMult"=rep(betaMu,2),"ncomponents"=ncomp, "betaTauMult"=diag(2),
                   "dirichParm"=dirichParm, "randPrecision"=randPrecision)
  
  initialValues = list("betaGender"=c(reg$coefficients[2]),"betaBy"=c(reg$coefficients[3]), "betaAge"=c(reg$coefficients[4]),
                       "errPrecision"=c(1), 
                       "Eta"=rep(1/ncomp, ncomp),
                       "S"=initS(nsubjects, ncomp),
                       "randmu_unordered"=cbind(randIntQuantiles, randSlopeQuantiles))
  
  stochasticNodes = c("betaGender","betaBy", "betaAge", "errPrecision",
                      "randmu", "randomComp", "Eta", "S")
  
  return(runJagsModel(datanodes, initialValues, stochasticNodes, jagsmodel, ncomp, niter, nthin, nburnin, nchains))
}

fitModel_errPrecision = function(niter=10000, nthin=50, nburnin=200, nchains=1, jagsmodel, ncomp=3, randmu_max, 
                                 randPrecision_max){
  
  dirichParm = rep(1, ncomp)
  reg=lm(weight~gender+by+age+time,data=ds)
  
  randPrecision = array(data=NA, dim = c(2,2,ncomp))
  randmu = matrix(data=NA, nrow=ncomp, ncol=2)
  for(i in 1:ncomp){
    randPrecision[,,i] = randPrecision_max[[i]]
    randmu[i,] = randmu_max[[i]]
  }
  
  datanodes = list("dsweight"=ds$weight,"dsby"=as.numeric(ds$by)-1, "dsage"=ds$age,
                   "dsgender"=as.numeric(ds$gender)-1,"dstime"=ds$time,
                   "nsubjects"=nsubjects,"nrep"=nrep, "gammaShapeRate"=gammaShapeRate, "betaMu"=betaMu,
                   "betaTau"=betaTau, "betaMuMult"=rep(betaMu,2),"ncomponents"=ncomp, "betaTauMult"=diag(2),
                   "dirichParm"=dirichParm, "randPrecision"=randPrecision, "randmu"=randmu)
  
  initialValues = list("betaGender"=c(reg$coefficients[2]),"betaBy"=c(reg$coefficients[3]), "betaAge"=c(reg$coefficients[4]),
                       "errPrecision"=c(1), 
                       "Eta"=rep(1/ncomp, ncomp),
                       "S"=initS(nsubjects, ncomp))
  
  stochasticNodes = c("betaGender","betaBy", "betaAge", "errPrecision",
                      "randomComp", "Eta", "S")
  
  return(runJagsModel(datanodes, initialValues, stochasticNodes, jagsmodel, ncomp, niter, nthin, nburnin, nchains))
}

fitModel_beta = function(niter=10000, nthin=50, nburnin=200, nchains=1, jagsmodel, ncomp=3, randmu_max, 
                                 randPrecision_max, errPrecision_max){
  
  dirichParm = rep(1, ncomp)
  reg=lm(weight~gender+by+age+time,data=ds)
  
  randPrecision = array(data=NA, dim = c(2,2,ncomp))
  randmu = matrix(data=NA, nrow=ncomp, ncol=2)
  for(i in 1:ncomp){
    randPrecision[,,i] = randPrecision_max[[i]]
    randmu[i,] = randmu_max[[i]]
  }
  
  datanodes = list("dsweight"=ds$weight,"dsby"=as.numeric(ds$by)-1, "dsage"=ds$age,
                   "dsgender"=as.numeric(ds$gender)-1,"dstime"=ds$time,
                   "nsubjects"=nsubjects,"nrep"=nrep, "gammaShapeRate"=gammaShapeRate, "betaMu"=betaMu,
                   "betaTau"=betaTau, "betaMuMult"=rep(betaMu,2),"ncomponents"=ncomp, "betaTauMult"=diag(2),
                   "dirichParm"=dirichParm, "randPrecision"=randPrecision, "randmu"=randmu, "errPrecision"=errPrecision_max)
  
  initialValues = list("betaGender"=c(reg$coefficients[2]),"betaBy"=c(reg$coefficients[3]), "betaAge"=c(reg$coefficients[4]),
                       "Eta"=rep(1/ncomp, ncomp),
                       "S"=initS(nsubjects, ncomp))
  
  stochasticNodes = c("betaGender","betaBy", "betaAge",
                      "randomComp", "Eta", "S")
  
  return(runJagsModel(datanodes, initialValues, stochasticNodes, jagsmodel, ncomp, niter, nthin, nburnin, nchains))
}

fitModel_eta = function(niter=10000, nthin=50, nburnin=200, nchains=1, jagsmodel, ncomp=3, randmu_max, 
                         randPrecision_max, errPrecision_max, betaBy_max, betaGender_max, betaAge_max){
  
  dirichParm = rep(1, ncomp)
  
  randPrecision = array(data=NA, dim = c(2,2,ncomp))
  randmu = matrix(data=NA, nrow=ncomp, ncol=2)
  for(i in 1:ncomp){
    randPrecision[,,i] = randPrecision_max[[i]]
    randmu[i,] = randmu_max[[i]]
  }
  
  datanodes = list("dsweight"=ds$weight,"dsby"=as.numeric(ds$by)-1, "dsage"=ds$age,
                   "dsgender"=as.numeric(ds$gender)-1,"dstime"=ds$time,
                   "nsubjects"=nsubjects,"nrep"=nrep, "gammaShapeRate"=gammaShapeRate, "betaMu"=betaMu,
                   "betaTau"=betaTau, "betaMuMult"=rep(betaMu,2),"ncomponents"=ncomp, "betaTauMult"=diag(2),
                   "dirichParm"=dirichParm, "randPrecision"=randPrecision, "randmu"=randmu, "errPrecision"=errPrecision_max,
                   "betaGender"=betaGender_max,"betaBy"=betaBy_max, "betaAge"=betaAge_max)
  
  initialValues = list("Eta"=rep(1/ncomp, ncomp), "S"=initS(nsubjects, ncomp))
  
  stochasticNodes = c("randomComp", "Eta", "S")
  
  return(runJagsModel(datanodes, initialValues, stochasticNodes, jagsmodel, ncomp, niter, nthin, nburnin, nchains))
}

