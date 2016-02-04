########## BAYESIAN PACKAGES ###########
install.packages('R2jags')
install.packages("mcmcplots")
install.packages("superdiag")
install.packages('mcmcse')
install.packages('ggmcmc')
library(R2jags)
library(mcmcplots)
library(superdiag)
library(mcmcse)
library(ggmcmc)

#make 4 groups
#weight is gender+by+quadtime+random intercept trimodal
#random intercept corresponds to diet---poor, ok, good
#time is 1 to 10
#males are 40 kg extra, 1997 is 30kg extra, quadtime coeff = 2, intercept is 50

time = 1:10
weightgen=function(gender, by){
  weight = rep(50, length(time));
  
  if(gender=="M"){
    weight = weight + 40
  }
  
  if(by=="97"){
    weight = weight + 30
  }
  
  weight = weight + time*10
  
  #random effect
  weight = weight + rnorm(1, 0, 10)
  
  weight = weight + rnorm(length(time), 0, 3)
}

#60 males, 60 females, 30 each from 96 and 97
sublist = matrix(nrow = 0, ncol = 5)
for(gender in c("M", "F")){
  for(by in c("96", "97")){
      for(k in 1:15){
        id = rep(paste(gender, by, k, sep=""), length(time))
        weight = weightgen(gender, by)
        subject = cbind(id, rep(gender, length(time)), rep(by, length(time)), 
                        time, weight)
        sublist=rbind(sublist, subject)
    }
  }
}

ds = data.frame(sublist)
colnames(ds) = c("id", "gender", "by", "time", "weight")
ds$time =  as.numeric(as.character(ds$time))
ds$weight = as.numeric(as.character(ds$weight))
ds = ds[order(ds$id, ds$time),]

########## BAYESIAN GLMM ANALYSIS ###########
nobs = length(ds$id)
nfolks = nobs/length(time)
nrep = length(time)
betaPrecision = 0.00001

dsweight = ds$weight
dstime = ds$time
dsgender = as.numeric(ds$gender)-1
dsby = as.numeric(ds$by)-1

model=function(){
  for(i in 0:(nfolks-1)){
    for(j in 1:nrep){
      dsweight[i*nrep+j]~dnorm(mu[i*nrep+j],errPrecision)
      mu[i*nrep+j]<-intercept+randomIntercept[i+1]+betaGender*dsgender[i*nrep+j]+betaBy*dsby[i*nrep+j]+betaTime*dstime[i*nrep+j]
    }
    randomIntercept[i+1]~dnorm(0,randPrecision)
  }
  
  randPrecision~dgamma(0.0005, 0.0005)
  errPrecision~dgamma(0.0005,0.0005)
  
  #Low precision is second parameter
  intercept~dnorm(0,betaPrecision)
  betaGender~dnorm(0,betaPrecision)
  betaBy~dnorm(0,betaPrecision)
  betaTime~dnorm(0,betaPrecision)
}

datanodes = c("dsweight","dsby","dsgender","dstime","nobs","nfolks", "nrep", "betaPrecision")
initialValues = list(list("intercept"=c(0),"betaGender"=c(0),"betaBy"=c(0),"betaTime"=c(0),"errPrecision"=c(1), "randPrecision"=c(1)))
stochasticNodes = c("intercept","betaGender","betaBy","betaTime","errPrecision", "randPrecision")

unload.module("glm")
fit = jags(data=datanodes, inits=initialValues, stochasticNodes, 
               n.chains=1, n.iter=100000,n.thin=5, n.burnin=30, 
               model.file=model, jags.module=NULL)
mcmcfit = as.mcmc(fit)

for(i in stochasticNodes){
  if(grepl("Preci+", i)){
    print(paste(i, sqrt(1/mean(mcmcfit[,i]))))
  }else{
    print(paste(i, mean(mcmcfit[,i])))
  }
}
