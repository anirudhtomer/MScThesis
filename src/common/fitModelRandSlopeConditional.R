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
wishartPriorDf = 4

fitModel = function(niter=10000, nthin=50, nburnin=200, nchains=1, jagsmodel, ncomp=3){
  
  dirichParm = rep(1, ncomp)
  randComp=extractRandomComp(viaReg = T)
  quantiles = seq(1/(ncomp+1),ncomp/(ncomp+1),length.out = ncomp)
  randIntQuantiles = quantile(randComp[,1], probs = quantiles)
  randSlopeQuantiles = quantile(randComp[,2], probs = quantiles)
  reg=lm(weight~gender+by+time,data=ds)
  
  datanodes = list("dsweight"=ds$weight,"dsby"=as.numeric(ds$by)-1,
                   "dsgender"=as.numeric(ds$gender)-1,"dstime"=ds$time,
                   "nsubjects"=nsubjects,"nrep"=nrep, "gammaShapeRate"=gammaShapeRate, "betaMu"=betaMu,
                   "betaTau"=betaTau, "betaMuMult"=rep(betaMu,2),"ncomponents"=ncomp, "betaTauMult"=diag(2),
                  "dirichParm"=dirichParm, "wishartPriorScale"=wishartPriorScale, "wishartPriorDf"=wishartPriorDf)
  
  initialValues = list("betaGender"=c(reg$coefficients[2]),"betaBy"=c(reg$coefficients[3]),
                       "errPrecision"=c(1), "precision1"=rep(1, ncomp), 
                       "precision2"=rep(1, ncomp),"rho"=rep(0, ncomp),
                       "Eta"=rep(1/ncomp, ncomp),
                       "S"=initS(nsubjects, ncomp),
                       "randmu"=cbind(randIntQuantiles, randSlopeQuantiles))
  #"randPrecision"=array(c(solve(var(splaash1))), c(2,2,ncomp)))
  
  stochasticNodes = c("betaGender","betaBy", "errPrecision", "randSigma","randPrecision",
                      "randmu", "randomComp", "Eta", "S", "precisionIntercept","precisionSlope", "rho")
  
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
  return(fit)
}