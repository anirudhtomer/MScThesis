install.packages('R2jags')
install.packages("moments")
install.packages('ggmcmc')
install.packages('MVN')
install.packages('beepr')
library(R2jags)
library(ggmcmc)
library(doParallel)
library(moments)
library(MASS)
library(MVN)
library(beepr)

registerDoParallel(cores = 8)

############ MIXTURE WITH INTERCEPT 1st way ###############
sigma <- matrix(c(16,11.9,11.9,9),2,2)

nobs = 90
nobsPerComp=c(floor(nobs/3),floor(nobs/3),floor(nobs/3))
first=mvrnorm(nobsPerComp[1], c(-30,-20.7), sigma)
second=mvrnorm(nobsPerComp[2], c(0,0), sigma)
third=mvrnorm(nobsPerComp[3], c(30,20.7), sigma)
sample = rbind(first, second, third)
nobs = nrow(sample)

qplot(x = sample[,1],y=sample[,2], color=
        unlist(lapply(nobsPerComp, function(i){rep(paste("C",i,sep=""), i)})))

ndim=2
model=function(){
  for(i in 1:nobs){
    sample[i,]~dmnorm(mu[S[i],1:ndim], omega[1:ndim, 1:ndim, S[i]])
    S[i]~dcat(eta[])
  }
  
  for(j in 1:ncomponents){
    mu[j, 1:ndim]~dmnorm(betaMuMult, betaTauMult)
    #mu[j,1]~dnorm(0,0.0001)
    #mu[j,2]~dnorm(0,0.0001)
   
    # omega[1:ndim, 1:ndim, j]~dwish(wishScale, 3)
    # sigma[1:ndim, 1:ndim, j]<-inverse(omega[1:ndim, 1:ndim, j])
    
    cor[j]~dunif(-1,1)
    precision1[j]~dgamma(0.00001, 0.00001)
    precision2[j]~dgamma(0.00001, 0.00001)

    sigma[1,1,j]<-1/precision1[j]
    sigma[2,2,j]<-1/precision2[j]
    sigma[1,2,j]<-cor[j] * sqrt(1/precision1[j]) * sqrt(1/precision2[j])
    sigma[2,1,j]<-cor[j] * sqrt(1/precision1[j]) * sqrt(1/precision2[j])

    omega[1:ndim,1:ndim,j]<-inverse(sigma[1:ndim, 1:ndim,j])
  }
  
  eta~ddirch(dirichParm[])
}

ncomponents=3
dirichParm = rep(1, ncomponents)

betaMuMult = c(0,0)
betaTauMult = diag(2)*0.0001
wishScale = diag(2)*0.0001
datanodes = c("sample","nobs","ncomponents", "dirichParm", "ndim", "betaMuMult", "betaTauMult", "wishScale")
tryMeanDim1=quantile(sample[,1], probs = seq(1/(ncomponents+1),ncomponents/(ncomponents+1), length.out = ncomponents))
tryMeanDim2=quantile(sample[,2], probs = seq(1/(ncomponents+1),ncomponents/(ncomponents+1), length.out = ncomponents))
initialValues = list(list("precision1"=rep(1, ncomponents), "precision2"=rep(1, ncomponents),
                          "eta"=rep(1/ncomponents, ncomponents), 
                          "mu"=matrix(c(tryMeanDim1, tryMeanDim2), nrow = ncomponents, ncol = ndim, byrow = F)))
params = c("sigma","omega","mu","eta", "S", "cor")

unload.module("glm")
nthin = 10
niter = 4000
nburnin = 0
fit = jags(data=datanodes, inits=initialValues, params, 
           n.chains=1, n.iter=niter,n.thin=nthin, n.burnin=0, 
           model.file=model, jags.module=NULL)

mcmcfit = as.mcmc(fit)
mcmcfit_partial = list(mcmcfit[[1]][((nburnin/nthin+1):(niter/nthin)),])
attributes(mcmcfit_partial) = attributes(mcmcfit)
mcmcLen = nrow(mcmcfit_partial[[1]])
attributes(mcmcfit_partial[[1]])$mcpar = c(1, mcmcLen, 1)

matrix(c(mean(mcmcfit[[1]][,"sigma[1,1,1]"]), mean(mcmcfit[[1]][,"sigma[1,2,1]"]),mean(mcmcfit[[1]][,"sigma[1,2,1]"]),mean(mcmcfit[[1]][,"sigma[2,2,1]"])), nrow=2, ncol=2)
matrix(c(mean(mcmcfit[[1]][,"sigma[1,1,2]"]), mean(mcmcfit[[1]][,"sigma[1,2,2]"]),mean(mcmcfit[[1]][,"sigma[1,2,2]"]),mean(mcmcfit[[1]][,"sigma[2,2,2]"])), nrow=2, ncol=2)
matrix(c(mean(mcmcfit[[1]][,"sigma[1,1,3]"]), mean(mcmcfit[[1]][,"sigma[1,2,3]"]),mean(mcmcfit[[1]][,"sigma[1,2,3]"]),mean(mcmcfit[[1]][,"sigma[2,2,3]"])), nrow=2, ncol=2)

matrix(c(mean(temp[[1]][,"sigma[1,1,1]"]), mean(temp[[1]][,"sigma[1,2,1]"]),mean(temp[[1]][,"sigma[1,2,1]"]),mean(temp[[1]][,"sigma[2,2,1]"])), nrow=2, ncol=2)
matrix(c(mean(temp[[1]][,"sigma[1,1,2]"]), mean(temp[[1]][,"sigma[1,2,2]"]),mean(temp[[1]][,"sigma[1,2,2]"]),mean(temp[[1]][,"sigma[2,2,2]"])), nrow=2, ncol=2)
matrix(c(mean(temp[[1]][,"sigma[1,1,3]"]), mean(temp[[1]][,"sigma[1,2,3]"]),mean(temp[[1]][,"sigma[1,2,3]"]),mean(temp[[1]][,"sigma[2,2,3]"])), nrow=2, ncol=2)

ggsobject = ggs(mcmcfit)
ggs_density(ggsobject, "mu")
ggs_density(ggsobject, "sigma")
  ggs_compare_partial(ggsobject,"mu")
ggs_running(ggsobject, "mu")

mcmcLen = (niter-nburnin)/nthin

multMean_DetCov=foreach(i=1:mcmcLen, .combine='rbind', .packages=c('MASS', 'MVN')) %do%{
  eta = getEta(mcmcfit, mcmcIterNum = i)
  randmu = getRandMu(mcmcfit, mcmcIterNum = i)
  sigmaEst = getSigma(mcmcfit, mcmcIterNum = i)

  postPred = matrix(nrow=0, ncol=2)
  for(j in 1:length(eta)){
    sampleSize = floor(1000*eta[j])
    if(sampleSize>0){
      postPred=rbind(postPred, mvrnorm(sampleSize, randmu[j,], sigmaEst[[j]]))
    }
  }
  
  return(c(apply(postPred, MARGIN=2, FUN = mean), det(cov(postPred))))
}
postPredDetCov=multMean_DetCov[,3]
sum(postPredDetCov>det(cov((sample))))

mahalaDist=foreach(i=1:mcmcLen, .combine='rbind', .packages=c('MASS', 'MVN')) %do%{
  eta = getEta(mcmcfit, mcmcIterNum = i)
  randmu = getRandMu(mcmcfit, mcmcIterNum = i)
  sigmaEst = getSigma(mcmcfit, mcmcIterNum = i)
  
  obsMahal = rep(0, nobs)
  for(j in 1:nobs){
    alloc = mcmcfit[[1]][i, paste("S[",j, "]", sep="")]
    for(m in 1:ncomponents){
      if(alloc!=m){
        obsMahal[j] = obsMahal[j] + mahalanobis(sample[j,], randmu[m,], sigmaEst)
      }
    }
    obsMahal[j] = obsMahal[j]/(ncomponents-1)
  }
  #obsDiff = c(mean(obsDiff[,1]), mean(obsDiff[,2]))
  
  return(mean(obsMahal))
  #newDiff = c(mean(newDiff[,1]), mean(newDiff[,2]))
  
  #newDiffMardia = attributes(mardiaTest(newDiff))
  #obsDiffMardia = attributes(mardiaTest(obsDiff))
  #return(c(newDiffMardia$g2p>obsDiffMardia$g2p,
  #        newDiffMardia$g1p>obsDiffMardia$g1p,
  #     attributes(hzTest(newDiff))$HZ>attributes(hzTest(obsDiff))$HZ))
}
beep(sound=8)
HPDinterval(mcmc(mahalaDist))
sum(mahalaDist[,1])/mcmcLen
sum(mahalaDist[,2])/mcmcLen
sum(mahalaDist[,3])/mcmcLen
HPDinterval(mcmc(mahalaDist))
plot(density(mahalaDist))

postPredictive=foreach(i=1:mcmcLen, .combine='rbind', .packages = 'moments') %dopar%{
  eta = getEta(mcmcfit, mcmcIterNum = i)
  randmu = getRandMu(mcmcfit, mcmcIterNum = i)
  precision = mcmcfit[[1]][i, "precision"]
  allocations = rmultinom(1, nobs, eta)
  
  postPred = numeric()
  for(j in 1:length(eta)){
    postPred = c(postPred, rnorm(n = allocations[j], sd = sqrt(1/precision), mean = randmu[j]))
  }
  postPred = postPred + intercept
  
  c(mean(postPred), var(postPred), skewness(postPred), kurtosis(postPred),
    IQR(postPred), max(postPred), min(postPred), moment(postPred, order = 6, central = T))
}

HPDinterval(mcmc(postPredictive[,1]))
mean(sample)
HPDinterval(mcmc(postPredictive[,2]))
var(sample)
HPDinterval(mcmc(postPredictive[,3]))
skewness(sample)
HPDinterval(mcmc(postPredictive[,4]))
kurtosis(sample)
HPDinterval(mcmc(postPredictive[,5]))
IQR(sample)
HPDinterval(mcmc(postPredictive[,6]))
max(sample)
HPDinterval(mcmc(postPredictive[,7]))
min(sample)
HPDinterval(mcmc(postPredictive[,8]))
moment(sample, order = 6, central = T)

densityplot = ggplot()+ aes(means) + geom_density()
densityplot + ylab(expression("p"[ Y ]*"(y)"))  + xlab("Random effect") + theme(axis.text=element_text(size=14),axis.title=element_text(size=18), plot.title=element_text(size=20))

getSigma = function(mcmcfit, mcmcIterNum){
  sigma = list(ncomponents)
  for(i in 1:ncomponents){
    sigma[[i]] = matrix(nrow = ndim, ncol = ndim)
    sigma[[i]][1,1] = mcmcfit[[1]][mcmcIterNum, paste("sigma[1,1,",i,"]",sep="")]
    sigma[[i]][1,2] = mcmcfit[[1]][mcmcIterNum, paste("sigma[1,2,",i,"]",sep="")]
    sigma[[i]][2,2] = mcmcfit[[1]][mcmcIterNum, paste("sigma[2,2,",i,"]",sep="")]
    sigma[[i]][2,1] = mcmcfit[[1]][mcmcIterNum, paste("sigma[2,1,",i,"]",sep="")]
  }
  return(sigma)
}

getRandMu=function(mcmcfit,  summaryFuncName="mean",  mcmcIterNum = NA){
  summaryFunc = get(summaryFuncName)
  randmu=matrix(nrow=ncomponents, ncol=2)
  
  for(i in 1:ncomponents){
    if(!is.numeric(mcmcIterNum)){
      randmu[i,] = c(summaryFunc(mcmcfit[[1]][,paste("mu[",i,",1]", sep="")]),
                     summaryFunc(mcmcfit[[1]][,paste("mu[",i,",2]", sep="")]))
    }else{
      randmu[i,] = c(mcmcfit[[1]][mcmcIterNum, paste("mu[",i,",1]", sep="")],
                     mcmcfit[[1]][mcmcIterNum, paste("mu[",i,",2]", sep="")])
    }
  }
  return(randmu)
}

getEta=function(mcmcfit, summaryFuncName="mean", mcmcIterNum = NA){
  summaryFunc = get(summaryFuncName)
  eta=numeric()
  
  for(i in 1:ncomponents){
    if(!is.numeric(mcmcIterNum)){
      eta[i] = summaryFunc(mcmcfit[[1]][,paste("eta[",i,"]", sep="")])
    }else{
      eta[i] = mcmcfit[[1]][mcmcIterNum, paste("eta[",i,"]", sep="")]
    }
  }
  
  return(eta)
}
