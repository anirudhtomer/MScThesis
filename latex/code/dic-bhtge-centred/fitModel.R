initS=function(nfolks, ncomponents){
  retval = numeric(nfolks)
  randomorder=sample(1:nfolks, nfolks)
  for(i in 1:ncomponents){
    retval[randomorder[((i-1)*(nfolks/ncomponents)+1):(i*(nfolks/ncomponents))]] = i;
  }
  return (retval)
}

fitModel = function(niter=10000, nthin=50, nburnin=200, nchains=1, jagsmodel, ncomp=3){
  datanodes = list("dsweight"=ds$weight,"dsby"=as.numeric(ds$by)-1,
                   "dsgender"=as.numeric(ds$gender)-1,"dstime"=ds$time,
                   "nsubjects"=nsubjects,"nrep"=nrep, 
                "varGammaParm"=0.0001, "betaTau"=0.0001, "betaMu"=0,
                "ncomponents"=ncomponents, "dirichParm"=rep(1, ncomp))
  initialValues = list("betaGender"=c(0),"betaBy"=c(0), "betaTime"=c(0),
                       "errPrecision"=c(1),
                       "randPrecision"=rep(1, 1),
                       "Eta"=rep(1/ncomp, ncomp),
                       "mue"=quantile(extractRandomComp(viaReg = T, nointercept = T), probs = seq(1/(ncomp+1),ncomp/(ncomp+1), length.out = ncomp)), 
                       "S"=initS(nsubjects, ncomp))
  stochasticNodes = c("betaGender","betaBy","betaTime", "errPrecision", 
                      "randPrecision", "randmu", "Eta", "S", "randomIntercept")
  
  if(ncomponents==1){
    #The following line removes the element Eta from initial values
    initialValues$Eta = NULL
    datanodes$dirichParm = NULL
    datanodes$S = initialValues$S
    initialValues$S=NULL
    datanodes$Eta=c(1)
    datanodes$ncomponents = NULL
    initialValues$mue=NULL
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