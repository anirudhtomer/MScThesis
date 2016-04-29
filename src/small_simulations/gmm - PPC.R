install.packages('R2jags')
install.packages("moments")
install.packages('ggmcmc')
library(R2jags)
library(ggmcmc)
library(doParallel)
library(moments)
library(beepr)

registerDoParallel(cores = 8)

############ MIXTURE WITH INTERCEPT 1st way ###############
nobs = 1000
first=rnorm(floor(nobs*0.6), -16, 3)
second=rnorm(floor(nobs*0.3), 0, 3)
third=rnorm(floor(nobs*0.1), 16, 3)
sample = c(first, second, third)
nobs = length(sample)
densityplot = ggplot()+ aes(sample) + geom_density()
densityplot + ylab(expression("p"[ Y ]*"(y)"))  + xlab("Random effect") + theme(axis.text=element_text(size=14),axis.title=element_text(size=18), plot.title=element_text(size=20))

model=function(){
  for(i in 1:nobs){
    sample[i]~dnorm(mu[S[i]], precision[S[i]])
    S[i]~dcat(eta[])
  }
  
  for(j in 1:ncomponents){
    mue[j]~dnorm(0, 0.0001) 
    precision[j]~dgamma(0.0005, 0.0005)
  }
  mu<-sort(mue[]-mean(mue[]))
  
  eta~ddirch(dirichParm[])
}

nobs = length(sample)
ncomponents=5
dirichParm = rep(1, ncomponents)

datanodes = c("sample","nobs","ncomponents", "dirichParm")
initialValues = list(list("precision"=rep(1,ncomponents), 
                          "eta"=rep(1/ncomponents, ncomponents),
                          "mue"=rep(0, ncomponents)))
params = c("precision","mu","eta", "S")

unload.module("glm")
nthin = 10
niter = 20000
nburnin = 8000
fit = jags(data=datanodes, inits=initialValues, params, 
           n.chains=1, n.iter=niter,n.thin=nthin, n.burnin=0, 
           model.file=model, jags.module=NULL)
mcmcfit = as.mcmc(fit)
mcmcfit_partial = list(mcmcfit[[1]][((nburnin/nthin+1):(niter/nthin)),])
attributes(mcmcfit_partial) = attributes(mcmcfit)
mcmcLen = nrow(mcmcfit_partial[[1]])
attributes(mcmcfit_partial[[1]])$mcpar = c(1, mcmcLen, 1)

ggsobject = ggs(mcmcfit_partial)
ggs_traceplot(ggsobject, "mu")
ggs_density(ggsobject, "eta")
ggs_density(ggsobject, "mu")
ggs_traceplot(ggsobject, "mu")
ggs_compare_partial(ggsobject,"mu")
ggs_running(ggsobject, "precision")

mis_count=foreach(i=1:mcmcLen, .combine='rbind') %do%{
  eta = getEta(mcmcfit, mcmcIterNum = i)
  randmu = getRandMu(mcmcfit, mcmcIterNum = i)
  randPrecision = getRandPrecision(mcmcfit, mcmcIterNum = i)
  allocations=getAllocations(mcmcfit_partial, i)
  
  count = 0
  for(j in 1:nobs){
    selfMahala = mahalanobis(sample[j], randmu[allocations[j]], 1/randPrecision[allocations[j]])
    for(m in 1:ncomponents){
      if(allocations[j]!=m){
        otherMahala = mahalanobis(sample[j], randmu[m], 1/randPrecision[m])
        if(selfMahala > otherMahala){
          count = count + 1
          break
        }
      }
    }
  }
  
  return(count)
}
beep(sound=8)
HPDinterval(mcmc(mis_count))
plot(density(mis_count))

mahalaDist=foreach(i=1:mcmcLen, .combine='rbind') %do%{
  eta = getEta(mcmcfit, mcmcIterNum = i)
  randmu = getRandMu(mcmcfit, mcmcIterNum = i)
  randPrecision = getRandPrecision(mcmcfit, mcmcIterNum = i)
  allocations=getAllocations(mcmcfit_partial, i)
  
  obsMahal = rep(0, nobs)
  for(j in 1:nobs){
    for(m in 1:ncomponents){
      if(allocations[j]!=m){
        obsMahal[j] = obsMahal[j] + 
          mahalanobis(sample[j], randmu[allocations[m]], sqrt(1/randPrecision[allocations[m]]))
      }
    }
    obsMahal[j] = obsMahal[j]/(ncomponents-1)
  }
  
  return(mean(obsMahal))
}
beep(sound=8)
HPDinterval(mcmc(mahalaDist))

postPredictive=foreach(i=1:mcmcLen, .combine='rbind', .packages = 'moments') %dopar%{
  eta = getEta(mcmcfit, mcmcIterNum = i)
  randmu = getRandMu(mcmcfit, mcmcIterNum = i)
  randPrecision = getRandPrecision(mcmcfit, mcmcIterNum = i)
  
  postPred = numeric()
  for(j in 1:length(eta)){
    postPred = c(postPred, rnorm(n = floor(nobs*eta[j]), sd = sqrt(1/randPrecision[j]), mean = randmu[j]))
  }
  
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

getRandPrecision=function(mcmcfit,  summaryFuncName="mean",  mcmcIterNum = NA){
  summaryFunc = get(summaryFuncName)
  randPrecision=numeric()
  
  for(i in 1:ncomponents){
    if(!is.numeric(mcmcIterNum)){
      randPrecision[i] = summaryFunc(mcmcfit[[1]][,paste("precision[",i,"]", sep="")])
    }else{
      randPrecision[i] = mcmcfit[[1]][mcmcIterNum, paste("precision[",i,"]", sep="")]
    }
  }
  randPrecision[randPrecision<1e-300]=1e-300
  return(randPrecision)
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

getAllocations=function(mcmcfit, mcmcIterNum){
  allocation = numeric(nobs)
  for(j in 1:nobs){
    allocation[j] = mcmcfit[[1]][mcmcIterNum, paste("S[",j,"]", sep="")]
  }
  return(allocation)
}
