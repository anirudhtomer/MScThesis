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
sigma <- matrix(c(16,3,3,9),2,2)

mat=rmultinom(1, 1000, c(1/4,1/4,1/4,1/4))
first=mvrnorm(mat[1,1], c(-108,0), sigma)
second=mvrnorm(mat[2,1], c(-7,0), sigma)
third=mvrnorm(mat[3,1], c(7,0), sigma)
fourth=mvrnorm(mat[4,1], c(18,0), sigma)
sample = rbind(first, second, third, fourth)

componentCol=c(rep("C1", mat[1,1]),rep("C2", mat[2,1]),
               rep("C3", mat[3,1]),rep("C4", mat[4,1]))

qplot(x = sample[,1],y=sample[,2], color=componentCol)

nobs = nrow(sample)
ncomponents=4
dirichParm = rep(1, ncomponents)

model=function(){
  for(i in 1:nobs){
    sample[i,]~dmnorm(mu[S[i],1:2], omega)
    S[i]~dcat(eta[])
  }
  
  for(j in 1:ncomponents){
    mu[j,1]~dnorm(0, 0.0001)
    mu[j,2]~dnorm(0, 0.0001)
  }
  
  sigma[1,1]<-var1
  sigma[2,2]<-var2
  sigma[1,2]<-cor * sqrt(var1) * sqrt(var2)
  sigma[2,1]<-cor * sqrt(var1) * sqrt(var2)
  
  var1<-1/precision1
  var2<-1/precision2
  cor~dunif(-1,1)
  
  omega<-inverse(sigma)
  
  precision1~dgamma(0.0005, 0.0005)
  precision2~dgamma(0.0005, 0.0005)
  eta~ddirch(dirichParm[])
}

datanodes = c("sample","nobs","ncomponents", "dirichParm")
initialValues = list(list("precision1"=c(1), "precision2"=c(1),
                          "eta"=rep(1/ncomponents, ncomponents)))
params = c("sigma","omega","mu","eta", "S")

unload.module("glm")
niter=10000
nthin=10
nburn=1000
fit = jags(data=datanodes, inits=initialValues, params, 
           n.chains=1, n.iter=niter,n.thin=nthin, n.burnin=nburn, 
           model.file=model, jags.module=NULL)
mcmcfit = as.mcmc(fit)
beep(sound=8)

ggsobject = ggs(mcmcfit)
ggs_density(ggsobject, "mu")
ggs_compare_partial(ggsobject,"mu")
ggs_running(ggsobject, "mu")

mcmcLen = (niter-nburn)/nthin

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
  sigma = matrix(nrow = 2, ncol = 2)
  sigma[1,1] = mcmcfit[[1]][mcmcIterNum, "sigma[1,1]"]
  sigma[1,2] = mcmcfit[[1]][mcmcIterNum, "sigma[1,2]"]
  sigma[2,1] = mcmcfit[[1]][mcmcIterNum, "sigma[2,1]"]
  sigma[2,2] = mcmcfit[[1]][mcmcIterNum, "sigma[2,2]"]
  
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
