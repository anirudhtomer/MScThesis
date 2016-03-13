initS=function(nfolks, ncomponents){
  retval = numeric(nfolks)
  randomorder=sample(1:nfolks, nfolks)
  for(i in 1:ncomponents){
    retval[randomorder[((i-1)*(nfolks/ncomponents)+1):(i*(nfolks/ncomponents))]] = i;
  }
  return (retval)
}

generateCovMatrixNonDiagonalIndices=function(nrep){
  indices = numeric()
  k=1;
  for(i in 1:nrep){
    for(j in 1:nrep){
      if(i!=j){
        indices[k] = i
        indices[k+1] = j
        k=k+2
      }
    }
  }
  return(indices)
}

getDirichParam = function(ncomp){
  rep(1, ncomp)
}

betaMu=0
betaTau=0.0001
gammaShapeRate = 0.0001

reg=lm(ds$weight~I(as.numeric(ds$gender)-1)+I(as.numeric(ds$by)-1)+ds$time)

initialBetaGender = reg$coefficients[2]
initialBetaBy = reg$coefficients[3]
initialBetaTime = reg$coefficients[4]

#use consecutive ones
fitModel = function(niter=10000, nthin=50, nburnin=200, nchains=1, jagsmodel, ncomp=3){
  
  datanodes = list("dsweight"=ds$weight,"dsby"= as.numeric(ds$by)-1,"dstime"=ds$time,
                   "dsgender"= as.numeric(ds$gender)-1,"nsubjects"=nrow(dswide), 
                   "nrep"=nrep, "betaMu"=betaMu, "betaTau"=betaTau,"ncomponents"=ncomp,
                   "gammaShape"=gammaShapeRate, "gammaRate"=gammaShapeRate,
                   "dirichParm"=getDirichParam(ncomp))
    
  initialValues = list("betaGender"=initialBetaGender,"betaBy"=initialBetaBy, "betaTime"=initialBetaTime,
                       "errPrecision"=c(1),
                       "randPrecision"=rep(1, 1),
                       "Eta"=rep(1/ncomp, ncomp),
                       "randmue"=quantile(extractRandomComp(viaReg = T), probs = seq(1/(ncomp+1),ncomp/(ncomp+1), length.out = ncomp)), 
                       "S"=initS(nsubjects, ncomp))
  
  stochasticNodes = c("XBeta", "randmu", "Eta", "S", 
                      "betaGender", "betaBy", "betaTime", "errPrecision",
                      "randomIntercept", "randPrecision")
  
  unload.module("glm")
  chainInits = list(initialValues)
  nchainsleft = nchains-1
  
  while(nchainsleft > 0){
    chainInits[[nchains - nchainsleft + 1]] = initialValues
    nchainsleft = nchainsleft - 1
  }
  
  fit = jags(data=datanodes, inits=chainInits, stochasticNodes,
             n.chains=nchains, n.iter=niter, n.thin=nthin, n.burnin=0,
             model.file=jagsmodel, jags.module=NULL)
  mcmcfit = as.mcmc(fit)
  mcmcfit[[1]] = mcmcfit[[1]][(nburnin/nthin+1):(niter/nthin),]
  
  return(list("fit" = fit, "mcmcfit" = mcmcfit))
}

fitModel_randmu = function(niter=10000, nthin=50, nburnin=200, nchains=1, jagsmodel, 
                           ncomp=3, rPrecision){
  
  datanodes = list("dsweight"=ds$weight,"dsby"= as.numeric(ds$by)-1,"dstime"=ds$time,
                   "dsgender"= as.numeric(ds$gender)-1,"nsubjects"=nrow(dswide), 
                   "nrep"=nrep, "betaMu"=betaMu, "betaTau"=betaTau,"ncomponents"=ncomp,
                   "gammaShape"=gammaShapeRate, "gammaRate"=gammaShapeRate,
                   "dirichParm"=getDirichParam(ncomp), "randPrecision"=rPrecision)
  
  initialValues = list("betaGender"=initialBetaGender,"betaBy"=initialBetaBy, "betaTime"=initialBetaTime,
                       "errPrecision"=c(1),
                       "Eta"=rep(1/ncomp, ncomp),
                       "randmue"=quantile(extractRandomComp(viaReg = T), probs = seq(1/(ncomp+1),ncomp/(ncomp+1), length.out = ncomp)), 
                       "S"=initS(nsubjects, ncomp))
  
  stochasticNodes = c("XBeta", "randmu", "Eta", "S", 
                      "betaGender", "betaBy", "betaTime", "errPrecision",
                      "randomIntercept")
  
  unload.module("glm")
  chainInits = list(initialValues)
  nchainsleft = nchains-1
  
  while(nchainsleft > 0){
    chainInits[[nchains - nchainsleft + 1]] = initialValues
    nchainsleft = nchainsleft - 1
  }
  
  fit = jags(data=datanodes, inits=chainInits, stochasticNodes,
             n.chains=nchains, n.iter=niter, n.thin=nthin, n.burnin=0,
             model.file=jagsmodel, jags.module=NULL)
  mcmcfit = as.mcmc(fit)
  mcmcfit[[1]] = mcmcfit[[1]][(nburnin/nthin+1):(niter/nthin),]
  
  return(list("fit" = fit, "mcmcfit" = mcmcfit))
}

fitModel_errPrecision = function(niter=10000, nthin=50, nburnin=200, nchains=1, jagsmodel, ncomp=3, 
                                 rMu,rPrecision){
  
  datanodes = list("dsweight"=ds$weight,"dsby"= as.numeric(ds$by)-1,"dstime"=ds$time,
                   "dsgender"= as.numeric(ds$gender)-1,"nsubjects"=nrow(dswide), 
                   "nrep"=nrep, "betaMu"=betaMu, "betaTau"=betaTau,"ncomponents"=ncomp,
                   "gammaShape"=gammaShapeRate, "gammaRate"=gammaShapeRate,
                   "dirichParm"=getDirichParam(ncomp), 
                   "randPrecision"=rPrecision, "randmu"=rMu)
  
  initialValues = list("betaGender"=initialBetaGender,"betaBy"=initialBetaBy, "betaTime"=initialBetaTime,
                       "errPrecision"=c(1),
                       "Eta"=rep(1/ncomp, ncomp),
                       "S"=initS(nsubjects, ncomp))
  
  stochasticNodes = c("XBeta", "Eta", "S", 
                      "betaGender", "betaBy", "betaTime", "errPrecision",
                      "randomIntercept")
  
  unload.module("glm")
  chainInits = list(initialValues)
  nchainsleft = nchains-1
  
  while(nchainsleft > 0){
    chainInits[[nchains - nchainsleft + 1]] = initialValues
    nchainsleft = nchainsleft - 1
  }
  
  fit = jags(data=datanodes, inits=chainInits, stochasticNodes,
             n.chains=nchains, n.iter=niter, n.thin=nthin, n.burnin=0,
             model.file=jagsmodel, jags.module=NULL)
  mcmcfit = as.mcmc(fit)
  mcmcfit[[1]] = mcmcfit[[1]][(nburnin/nthin+1):(niter/nthin),]
  
  return(list("fit" = fit, "mcmcfit" = mcmcfit))
}

fitModel_beta = function(niter=10000, nthin=50, nburnin=200, nchains=1, jagsmodel, ncomp=3, 
                                 rMu,rPrecision, ePrecision){
  
  datanodes = list("dsweight"=ds$weight,"dsby"= as.numeric(ds$by)-1,"dstime"=ds$time,
                   "dsgender"= as.numeric(ds$gender)-1,"nsubjects"=nrow(dswide), 
                   "nrep"=nrep, "betaMu"=betaMu, "betaTau"=betaTau,"ncomponents"=ncomp,
                   "dirichParm"=getDirichParam(ncomp),"errPrecision"=ePrecision,
                   "randPrecision"=rPrecision, "randmu"=rMu)
  
  initialValues = list("betaGender"=initialBetaGender,"betaBy"=initialBetaBy, "betaTime"=initialBetaTime,
                       "Eta"=rep(1/ncomp, ncomp),
                       "S"=initS(nsubjects, ncomp))
  
  stochasticNodes = c("XBeta", "Eta", "S", 
                      "betaGender", "betaBy", "betaTime",
                      "randomIntercept")
  
  unload.module("glm")
  chainInits = list(initialValues)
  nchainsleft = nchains-1
  
  while(nchainsleft > 0){
    chainInits[[nchains - nchainsleft + 1]] = initialValues
    nchainsleft = nchainsleft - 1
  }
  
  fit = jags(data=datanodes, inits=chainInits, stochasticNodes,
             n.chains=nchains, n.iter=niter, n.thin=nthin, n.burnin=0,
             model.file=jagsmodel, jags.module=NULL)
  mcmcfit = as.mcmc(fit)
  mcmcfit[[1]] = mcmcfit[[1]][(nburnin/nthin+1):(niter/nthin),]
  
  return(list("fit" = fit, "mcmcfit" = mcmcfit))
}

fitModel_eta = function(niter=10000, nthin=50, nburnin=200, nchains=1, jagsmodel, ncomp=3, 
                         rMu,rPrecision, ePrecision, bGender, bBy, bTime){
  
  datanodes = list("dsweight"=ds$weight,"dsby"= as.numeric(ds$by)-1,"dstime"=ds$time,
                   "dsgender"= as.numeric(ds$gender)-1,"nsubjects"=nrow(dswide), 
                   "nrep"=nrep, "ncomponents"=ncomp,
                   "dirichParm"=getDirichParam(ncomp),"errPrecision"=ePrecision,
                   "randPrecision"=rPrecision, "randmu"=rMu, "betaGender"=bGender,
                   "betaBy"=bBy, "betaTime"=bTime)
  
  initialValues = list("Eta"=rep(1/ncomp, ncomp),
                       "S"=initS(nsubjects, ncomp))
  
  stochasticNodes = c("XBeta", "Eta", "S",  "randomIntercept")
  
  unload.module("glm")
  chainInits = list(initialValues)
  nchainsleft = nchains-1
  
  while(nchainsleft > 0){
    chainInits[[nchains - nchainsleft + 1]] = initialValues
    nchainsleft = nchainsleft - 1
  }
  
  fit = jags(data=datanodes, inits=chainInits, stochasticNodes,
             n.chains=nchains, n.iter=niter, n.thin=nthin, n.burnin=0,
             model.file=jagsmodel, jags.module=NULL)
  mcmcfit = as.mcmc(fit)
  mcmcfit[[1]] = mcmcfit[[1]][(nburnin/nthin+1):(niter/nthin),]
  
  return(list("fit" = fit, "mcmcfit" = mcmcfit))
}

