install.packages('R2jags')
install.packages('ggmcmc')
library(R2jags)
library(ggmcmc)

source("DIC_functions.R")

mat=rmultinom(1, 1000, c(1/3,1/3,1/3))
first=rnorm(mat[1,1], -6, 3)
second=rnorm(mat[2,1], 0, 3)
third=rnorm(mat[3,1], 6, 3)
sample = c(first, second, third)

densityplot = ggplot()+ aes(sample) + geom_density()
densityplot + ylab(expression("p"[ Y ]*"(y)"))  + xlab("Y") + theme(axis.text=element_text(size=14),axis.title=element_text(size=18), plot.title=element_text(size=20))

nobs = length(sample)
ncomponents=3

model=function(){
  for(i in 1:nobs){
    sample[i]~dnorm(mu[S[i]], precision)
    S[i]~dcat(eta[])
    
    for(m in 1:ncomponents){
      obsLikelihood[i,m]<-eta[m]*pow(precision/(2*pi),0.5)*exp(-0.5 * pow(sample[i]-mu[m],2) * precision)
      condLikelihood[i,m]<-log(obsLikelihood[i,m]) * obsLikelihood[i,m]
    }
    
    obsL[i]<-sum(obsLikelihood[i,])
    condLl[i]<-sum(condLikelihood[i,]/obsL[i])
    obsLl[i]<-log(obsL[i])
    regLl[i]<-log(pow(precision/(2*pi),0.5)*exp(-0.5 * pow(sample[i]-mu[S[i]],2) * precision))
  }
  
  for(j in 1:ncomponents){
    mue[j]~dnorm(0, 0.0001)  
  }
  
  mu<-sort(mue)
  
  condDeviance <- sum(condLl) * (-2)
  obsDeviance <- sum(obsLl) * (-2)
  regDeviance <- sum(regLl) * (-2)
  
  precision~dgamma(0.0005, 0.0005)
  eta~ddirch(dirichParm[])
}

datanodes = list("sample"=sample,"nobs"=nobs,"ncomponents"=ncomponents, 
                 "dirichParm"=rep(1, ncomponents), "pi"=pi)
initialValues = list(list("precision"=c(1),
                          "eta"=rep(1/ncomponents, ncomponents),
                          "mue"=quantile(sample, probs = seq(1/(ncomponents+1),ncomponents/(ncomponents+1), length.out = ncomponents))))
params = c("precision","mu","eta", "obsDeviance","obsL", "regDeviance", "condDeviance", "S")

unload.module("glm")
niter = 10000
nthin = 50
nburnin = 2000
fit = jags(data=datanodes, inits=initialValues, params, 
           n.chains=1, n.iter=niter,n.thin=nthin, n.burnin=nburnin, 
           model.file=model, jags.module=NULL)
mcmcfit = as.mcmc(fit)

fit
calculateDIC7(mcmcfit, "regDeviance")

saveResult(ncomponents,
           calculateObsDIC1(mcmcfit, "obsDeviance")$dic,
           calculateObsDIC2(mcmcfit, "obsDeviance")$dic,
           calculateObsDIC3(mcmcfit, "obsDeviance","obsL")$dic,
           calculateDIC7(mcmcfit, "condDeviance")$dic)

ggsobject = ggs(mcmcfit)
ggs_density(ggsobject, "mu")
ggs_running(ggsobject, "Devi")
ggs_running(ggsobject, "mu")
ggs_compare_partial(ggsobject, "Devi")