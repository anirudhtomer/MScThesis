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


#use consecutive ones
fitModel = function(niter=10000, nthin=50, nburnin=200, nchains=1, jagsmodel, ncomp=3){
  reg=lm(ds$weight~I(as.numeric(ds$gender)-1)+I(as.numeric(ds$by)-1)+ds$time)
  
  initialBetaGender = reg$coefficients[2]
  initialBetaBy = reg$coefficients[3]
  initialBetaTime = reg$coefficients[4]
  
  datanodes = list("weight"=dswide_y,"dsby"= as.numeric(dswide$by)-1,"dstime"=time,
                   "dsgender"= as.numeric(dswide$gender)-1,"nsubjects"=nrow(dswide), 
                   "nrep"=nrep, "betaMu"=0, "betaTau"=0.0001,"ncomponents"=ncomp,
                   "varGammaParm"=0.0001,"dirichParm"=rep(1, ncomp),
                   "betaRandMu"=rep(0, ncomp), "betaRandMuTau"=diag(ncomp),
                   "covMatNonDiagIndices"=generateCovMatrixNonDiagonalIndices(nrep))
    
  initialValues = list("betaGender"=c(initialBetaGender),"betaBy"=c(initialBetaBy), 
                       "betaTime"=c(initialBetaTime),
                       "errPrecision"=c(1),"correlation"=rep(0, ncomponents),
                       "randPrecision"=rep(1, ncomponents),
                       "Eta"=rep(1/ncomp, ncomp),
                       "randmu"=c(130,145, 165),
                       #"randmu"=quantile(extractRandomComp(viaReg = T), probs = seq(1/(ncomp+1),ncomp/(ncomp+1), length.out = ncomp)), 
                       "S"=initS(nsubjects, ncomp))
  stochasticNodes = c("XBeta","sigma", "omega", "randmu", "Eta", "S", 
                      "betaBy","betaTime", "betaGender","randPrecision", "errPrecision", "correlation")
  
  if(ncomponents==1){
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
  return(fit)
}