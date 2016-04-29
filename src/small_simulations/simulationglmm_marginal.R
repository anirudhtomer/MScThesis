########## BAYESIAN PACKAGES ###########
install.packages('R2jags')
install.packages('ggmcmc')
library(R2jags)
library(ggmcmc)
library(reshape2)

#make 4 groups
#weight is gender+by+quadtime+random intercept trimodal
#random intercept corresponds to diet---poor, ok, good
#time is 1 to 10
#males are 40 kg extra, 1997 is 30kg extra, quadtime coeff = 2, intercept is 50

time = 1:5
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

dswide=data.frame(dcast(ds, id+gender+by ~ time, value.var = "weight"))
dswide_y=as.matrix(dswide[,c(-1,-2,-3)])
########## BAYESIAN GLMM ANALYSIS ###########
nobs = nrow(ds)

model=function(){
  for(i in 1:nfolks){
    weight[i,(1:nrep)]~dmnorm(fixedEff[i,(1:nrep)], inverse(sigma[,]))
    for(j in 1:nrep){
      fixedEff[i,j]<-intercept+betaGender*dsgender[i]+
        betaBy*dsby[i]+betaTime*dstime[j]
    }
  }
  
  for(m in 1:(nrep*(nrep-1))){
    sigma[covMatNonDiagIndices[m*2-1],covMatNonDiagIndices[m*2]] = 1/randPrecision
  }
  
  for(n in 1:nrep){
    sigma[n,n]<-1/randPrecision + 1/errPrecision
  }
  
  randPrecision~dgamma(0.0005, 0.0005)
  errPrecision~dgamma(0.0005,0.0005)
  
  intercept~dnorm(0,betaPrecision)
  betaGender~dnorm(0,betaPrecision)
  betaBy~dnorm(0,betaPrecision)
  betaTime~dnorm(0,betaPrecision)
}

generateCovMatrixNonDiagonalIndices=function(ncomp){
  indices = numeric()
  k=1;
  for(i in 1:ncomp){
    for(j in 1:ncomp){
      if(i!=j){
        indices[k] = i
        indices[k+1] = j
        k=k+2
      }
    }
  }
  return(indices)
}

#use consecutive ones
covMatrixIndices = generateCovMatrixNonDiagonalIndices(length(time))

datanodes = list("weight"=dswide_y,"dsby"= as.numeric(dswide$by)-1,"dstime"=time,
                 "dsgender"= as.numeric(dswide$gender)-1,"nfolks"=nrow(dswide), 
                 "nrep"=length(time), "betaPrecision"=0.0001, "covMatNonDiagIndices"=covMatrixIndices)
initialValues = list(list("intercept"=c(0),"betaGender"=c(0),"betaBy"=c(0),
                          "betaTime"=c(0),"errPrecision"=c(1), 
                          "randPrecision"=c(1)))
stochasticNodes = c("intercept","betaGender","betaBy","betaTime",
                    "randPrecision", "errPrecision", "sigma")
unload.module("glm")
fit = jags(data=datanodes, inits=initialValues, stochasticNodes, 
               n.chains=1, n.iter=20000,n.thin=50, n.burnin=1000, 
               model.file=model, jags.module=NULL)
mcmcfit = as.mcmc(fit)
ggsobj = ggs(mcmcfit)

ggs_density(ggsobj, "temp")
ggs_running(ggsobj, "Prec")
ggs_autocorrelation(ggsobj,"Prec")
ggs_compare_partial(ggsobj, "Prec")
