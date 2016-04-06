initS=function(nfolks, ncomponents){
  retval = numeric(nfolks)
  randomorder=sample(1:nfolks, nfolks)
  for(i in 1:ncomponents){
    retval[randomorder[((i-1)*(nfolks/ncomponents)+1):(i*(nfolks/ncomponents))]] = i;
  }
  return (retval)
}

fitModel = function(niter=10000, nthin=50, nburnin=200, nchains=1, jagsmodel, ncomp=3){

  randComp=extractRandomComp(viaReg = T)
  quantiles = seq(1/(ncomponents+1),ncomponents/(ncomponents+1),length.out = ncomponents)
  randIntQuantiles = quantile(randComp[,1], probs = quantiles)
  randSlopeQuantiles = quantile(randComp[,2], probs = quantiles)
  reg=lm(weight~gender+by+time,data=ds)
  
  datanodes = list("dsweight"=ds$weight,"dsby"=as.numeric(ds$by)-1,
                   "dsgender"=as.numeric(ds$gender)-1,"dstime"=ds$time,
                   "nsubjects"=nsubjects,"nrep"=nrep, "varGammaParm"=0.0001,
                   "betaTau"=0.0001, "betaMu"=0,"ncomponents"=ncomp, 
                  "dirichParm"=rep(1, ncomp), "wishartParm"=diag(2)*9)
  
  initialValues = list("betaGender"=c(reg$coefficients[2]),"betaBy"=c(reg$coefficients[3]),
                       "errPrecision"=c(1), "precision1"=rep(1, ncomponents), 
                       "precision2"=rep(1, ncomponents),"rho"=rep(0, ncomponents),
                       "Eta"=rep(1/ncomponents, ncomponents),
                       "S"=initS(nsubjects, ncomponents),
                       "randmu"=cbind(randIntQuantiles, randSlopeQuantiles))
                       #"randPrecision"=solve(matrix(c(7.74,3.18,3.18,16.75), nrow=2, ncol=2)))
  
  stochasticNodes = c("betaGender","betaBy", "errPrecision", "randSigma","randPrecision",
                      "randmu", "randomComp", "Eta", "S", "precision1","precision2", "rho")
  
  dsweight = ds$weight
  dstime = ds$time
  dsgender = as.numeric(ds$gender)-1
  dsby = as.numeric(ds$by)-1
  
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