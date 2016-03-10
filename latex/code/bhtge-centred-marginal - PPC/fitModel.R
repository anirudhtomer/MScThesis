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
  
  datanodes = list("weight"=dswide_y,"dsby"= as.numeric(dswide$by)-1,"dstime"=time,
                   "dsgender"= as.numeric(dswide$gender)-1,"nsubjects"=nrow(dswide), 
                   "nrep"=nrep, "betaMu"=0, "betaTau"=0.0001,"ncomponents"=ncomp,
                   "varGammaParm"=0.0001,"dirichParm"=rep(1, ncomp),
                   "covMatNonDiagIndices"=generateCovMatrixNonDiagonalIndices(nrep))
    
  initialValues = list("betaGender"=c(0),"betaBy"=c(0), "betaTime"=c(0),
                       "errPrecision"=c(1),
                       "randPrecision"=rep(1, 1),
                       "Eta"=rep(1/ncomp, ncomp),
                       "randmu"=quantile(extractRandomComp(viaReg = T), probs = seq(1/(ncomp+1),ncomp/(ncomp+1), length.out = ncomp)), 
                       "S"=initS(nsubjects, ncomp))
  stochasticNodes = c("XBeta","sigma", "omega", "randmu", "Eta", "S")
  
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