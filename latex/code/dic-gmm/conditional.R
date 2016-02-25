install.packages("doParallel")
library(doParallel)

############## Conditional deviance again ############
dev = 0
mcmcLen = ((niter-nburnin)/nthin)
for(i in 1:mcmcLen){
  precision = mcmcfit[[1]][i, "precision"]
  for(j in 1:nobs){
    clusterNum = round(mcmcfit[[1]][i,paste("S[",j,"]", sep="")])
    dev = dev - 2*log(mcmcfit[[1]][i,paste("eta[",clusterNum,"]",sep="")])
    dev = dev + log(2*pi) - log(precision)
    dev = dev + precision * ((sample[j]-mcmcfit[[1]][i,paste("mu[",clusterNum,"]",sep="")])^2)
  }
}
dev=dev/mcmcLen

allocation = matrix(nrow=mcmcLen, ncol=length(sample))
for(i in 1:mcmcLen){
  for(j in 1:length(sample)){
    allocation[i,j] = round(mcmcfit[[1]][i,paste("S[",j,"]", sep="")])
  }
}

########## Conditional model ########
conditionalModel=function(){
  for(i in 1:nobs){
    sample[i]~dnorm(mu[S[i]], precision)
    S[i]~dcat(eta[])
  }    
  
  for(j in 1:ncomponents){
    mue[j]~dnorm(0, 0.0001)  
  }
  
  mu<-sort(mue)
  
  precision~dgamma(0.0005, 0.0005)
  eta~ddirch(dirichParm[])
}

registerDoParallel(cores=8)

res=foreach(i=1:mcmcLen) %dopar% {
  library(R2jags)
  library(ggmcmc)
  datanodes = list("sample"=sample,"nobs"=nobs,"ncomponents"=ncomponents, 
                   "dirichParm"=rep(1, ncomponents), "S" = allocation[1,])
  initialValues = list(list("precision"=c(1),
                            "eta"=rep(1/ncomponents, ncomponents),
                            "mue"=quantile(sample, probs = seq(1/(ncomponents+1),ncomponents/(ncomponents+1), length.out = ncomponents))))
  params = c("precision","mu","eta")
  
  unload.module("glm")
  niter = 10000
  nthin = 50
  nburnin = 2000
  fit = jags(data=datanodes, inits=initialValues, params, 
             n.chains=1, n.iter=niter,n.thin=nthin, n.burnin=nburnin, 
             model.file=conditionalModel, jags.module=NULL)
  mcmcfit = as.mcmc(fit)
}


## Find second part of the likelihood ###
dthetabar  = 0
for(i in 1:length(res)){
  tempmcmcfit = res[[i]][[1]]
  tempallocation = allocation[i,]
  precision = mean(tempmcmcfit[, "precision"])
  for(j in 1:nobs){
    clusterNum = tempallocation[j]
    dthetabar = dthetabar - 2*log(mean(tempmcmcfit[,paste("eta[",clusterNum,"]",sep="")]))
    dthetabar = dthetabar + log(2*pi) - log(precision)
    dthetabar = dthetabar + precision * ((sample[j]-mean(tempmcmcfit[,paste("mu[",clusterNum,"]",sep="")]))^2)
  }
}
dthetabar=dthetabar/mcmcLen
  
