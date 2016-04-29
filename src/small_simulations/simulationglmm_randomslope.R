########## BAYESIAN PACKAGES ###########
install.packages('R2jags')
install.packages('ggmcmc')
library(R2jags)
library(ggmcmc)
library(MASS)

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
  randeff=mvrnorm(1,c(0,0),matrix(c(100,120,120,400),2,2))
  weight = weight + randeff[1]
  weight = weight + time * randeff[2]
  
  weight = weight + rnorm(length(time), 0, 3)
}

#60 males, 60 females, 30 each from 96 and 97
sublist = matrix(nrow = 0, ncol = 5)
for(gender in c("M", "F")){
  for(by in c("96", "97")){
      for(k in 1:50){
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

randMuMean = c(0,0)
wishartParm=matrix(c(1,0,0,1),2,2)

model=function(){
  for(i in 0:(nfolks-1)){
    for(j in 1:nrep){
      dsweight[i*nrep+j]~dnorm(mu[i*nrep+j],errPrecision)
      mu[i*nrep+j]<-randomMu[i+1,1]+
        betaGender*dsgender[i*nrep+j]+betaBy*dsby[i*nrep+j]+
        randomMu[i+1,2]*dstime[i*nrep+j]
    }
    randomMu[i+1,(1:2)]~dmnorm(randMu[], sigmainv)
  }
  
  sigmainv~dwish(wishartParm[,],2)
  sigma<-inverse(sigmainv)
  errPrecision~dgamma(0.0005,0.0005)
  
  #Low precision is second parameter
  randMu[1]~dnorm(0,betaPrecision)
  betaGender~dnorm(0,betaPrecision)
  betaBy~dnorm(0,betaPrecision)
  randMu[2]~dnorm(0,betaPrecision)
}

datanodes = c("dsweight","dsby","dsgender","dstime","nfolks", 
              "nrep", "betaPrecision", "wishartParm")
initialValues = list(list("randMu"=c(0,0),"betaGender"=c(0),
                          "betaBy"=c(0),
                          "errPrecision"=c(1)))
stochasticNodes = c("betaGender","betaBy","randMu",
                    "errPrecision", "sigma", "sigmainv")

unload.module("glm")
fit = jags(data=datanodes, inits=initialValues, stochasticNodes, 
               n.chains=1, n.iter=100000,n.thin=100, n.burnin=3000, 
               model.file=model, jags.module=NULL)
mcmcfit = as.mcmc(fit)
ggsobj = ggs(mcmcfit)

ggs_density(ggsobj, "sigma")
ggs_pairs(ggsobj, "Precision")
