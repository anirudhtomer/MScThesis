initS=function(nfolks, ncomponents){
  retval = numeric(nfolks)
  randomorder=sample(1:nfolks, nfolks)
  for(i in 1:ncomponents){
    retval[randomorder[((i-1)*(nfolks/ncomponents)+1):(i*(nfolks/ncomponents))]] = i;
  }
  return (retval)
}

fitModel = function(niter=10000, nthin=50, nburnin=200, nchains=1, jagsmodel){
  betaTau = 0.0001
  betaMu = 0
  varGammaParm = 0.0001
  dirichParm = rep(1, ncomponents)
    
  datanodes = c("dsweight","dsby","dsgender","dstime","nsubjects","nrep", 
                "varGammaParm","betaTau", "betaMu","ncomponents", "dirichParm")
  initialValues = list("betaGender"=c(0),"betaBy"=c(0), "betaTime"=c(0),
                       "errPrecision"=c(1),"randPrecision"=rep(1, 1),
                       "Eta"=rep(1/ncomponents, ncomponents),
                       "mue"=quantile(extractRandomComp(viaReg = T), probs = seq(1/(ncomponents+1),ncomponents/(ncomponents+1), length.out = ncomponents)), 
                       "S"=initS(nsubjects, ncomponents))
  stochasticNodes = c("betaGender","betaBy",
                      "betaTime","errPrecision", "randPrecision", 
                      "randmu", "Eta")
  
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