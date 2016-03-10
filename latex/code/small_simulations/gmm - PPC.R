install.packages('R2jags')
install.packages("moments")
install.packages('ggmcmc')
library(R2jags)
library(ggmcmc)
library(doParallel)
library(moments)

registerDoParallel(cores = 8)

############ MIXTURE WITH INTERCEPT 1st way ###############
mat=rmultinom(1, 1000, c(1/3,1/3,1/3))
first=rnorm(mat[1,1], -16, 3)
second=rnorm(mat[2,1], 0, 3)
third=rnorm(mat[3,1], 16, 3)
sample = c(first, second, third)
densityplot = ggplot()+ aes(sample) + geom_density()
densityplot + ylab(expression("p"[ Y ]*"(y)"))  + xlab("Random effect") + theme(axis.text=element_text(size=14),axis.title=element_text(size=18), plot.title=element_text(size=20))

intercept = 20
sample = sample + intercept

nobs = length(sample)
ncomponents=2
dirichParm = rep(1, ncomponents)

model=function(){
  for(i in 1:nobs){
    sample[i]~dnorm(intercept + mu[S[i]], precision)
    S[i]~dcat(eta[])
  }
  
  for(j in 1:ncomponents){
    mue[j]~dnorm(0, 0.0001)  
  }
  mu<-sort(mue[]-mean(mue[]))
  
  precision~dgamma(0.0005, 0.0005)
  intercept~dnorm(0, 0.0001)
  eta~ddirch(dirichParm[])
}

datanodes = c("sample","nobs","ncomponents", "dirichParm")
initialValues = list(list("precision"=c(1), 
                          "eta"=rep(1/ncomponents, ncomponents),
                          "mue"=rep(0, ncomponents),
                          "intercept"=c(0)))
params = c("precision","mu","eta","intercept", "S")

unload.module("glm")
fit = jags(data=datanodes, inits=initialValues, params, 
           n.chains=1, n.iter=10000,n.thin=10, n.burnin=3000, 
           model.file=model, jags.module=NULL)
mcmcfit = as.mcmc(fit)

ggsobject = ggs(mcmcfit)
ggs_density(ggsobject, "mu")
ggs_compare_partial(ggsobject,"mu")
ggs_running(ggsobject, "mu")

mcmcLen = (3000-300)/5

mapeDiff=foreach(i=1:mcmcLen, .combine='c') %dopar%{
  eta = getEta(mcmcfit, mcmcIterNum = i)
  randmu = getRandMu(mcmcfit, mcmcIterNum = i)
  precision = mcmcfit[[1]][i, "precision"]
  
  mapeObs = 0
  mapeNew = 0
  for(j in 1:nobs){
    alloc = mcmcfit[[1]][i, paste("S[",j, "]", sep="")]
    mapeObs = mapeObs + abs(sample[j]-intercept-randmu[alloc])
    mapeNew = mapeNew + abs(rnorm(1, randmu[alloc], sqrt(1/precision))-randmu[alloc])
  }
  mapeObs = mapeObs/nobs
  mapeNew = mapeNew/nobs
  
  mapeNew-mapeObs
}
HPDinterval(mcmc(mapeDiff))
plot(density(mapeDiff))

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
  
  c(mean(postPred), var(postPred), skewness(postPred),kurtosis(postPred),
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

getRandMu=function(mcmcfit,  summaryFuncName="mean",  mcmcIterNum = NA){
  summaryFunc = get(summaryFuncName)
  randmu=numeric()
  
  for(i in 1:ncomponents){
    if(!is.numeric(mcmcIterNum)){
      randmu[i] = summaryFunc(mcmcfit[[1]][,paste("mu[",i,"]", sep="")])
    }else{
      randmu[i] = mcmcfit[[1]][mcmcIterNum, paste("mu[",i,"]", sep="")]
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
