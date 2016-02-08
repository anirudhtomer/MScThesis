install.packages('R2jags')
install.packages('ggmcmc')
library(R2jags)
library(ggmcmc)

############ MIXTURE WITH INTERCEPT 1st way ###############
mat=rmultinom(1, 1000, c(1/3,1/3,1/3))
first=rnorm(mat[1,1], -16, 3)
second=rnorm(mat[2,1], 0, 3)
third=rnorm(mat[3,1], 16, 3)
sample = c(first, second, third)

densityplot = ggplot()+ aes(sample) + geom_density()
densityplot + ylab(expression("p"[ Y ]*"(y)"))  + xlab("Y") + theme(axis.text=element_text(size=14),axis.title=element_text(size=18), plot.title=element_text(size=20))

intercept = 20
sample = sample + intercept

nobs = length(sample)
ncomponents=3
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
params = c("precision","mu","eta","intercept")

unload.module("glm")
fit = jags(data=datanodes, inits=initialValues, params, 
           n.chains=1, n.iter=3000,n.thin=5, n.burnin=30, 
           model.file=model, jags.module=NULL)
mcmcfit = as.mcmc(fit)

ggsobject = ggs(mcmcfit)
ggs_density(ggsobject, "mu")

############ MIXTURE WITH INTERCEPT 2nd WAY ###############
mat=rmultinom(1, 1000, c(1/3,1/3,1/3))
first=rnorm(mat[1,1], -16, 3)
second=rnorm(mat[2,1], 0, 3)
third=rnorm(mat[3,1], 16, 3)
sample = c(first, second, third)

densityplot = ggplot()+ aes(sample) + geom_density()
densityplot + ylab(expression("p"[ Y ]*"(y)"))  + xlab("Y") + theme(axis.text=element_text(size=14),axis.title=element_text(size=18), plot.title=element_text(size=20))

intercept = 20
sample = sample + intercept

nobs = length(sample)
ncomponents=3
dirichParm = rep(1, ncomponents)

model=function(){
  for(i in 1:nobs){
    sample[i]~dnorm(intercept + mu[S[i]], precision)
    S[i]~dcat(eta[])
  }
  
  for(j in 1:ncomponents){
    mu[j]~dnorm(0, 0.0001)  
  }
  
  z~dbern(constraint)
  constraint <- step(abs(mu[1] + mu[2] + mu[3]) - 0.05)
  
  precision~dgamma(0.0005, 0.0005)
  intercept~dnorm(0, 0.0001)
  eta~ddirch(dirichParm[])
}
z = 0

datanodes = c("sample","nobs","ncomponents", "dirichParm", "z")
initialValues = list(list("precision"=c(1), 
                          "eta"=rep(1/ncomponents, ncomponents),
                          "mu"=rep(0, ncomponents),
                          "intercept"=c(0)))
params = c("precision","mu","eta","intercept")

unload.module("glm")
fit = jags(data=datanodes, inits=initialValues, params, 
           n.chains=1, n.iter=3000,n.thin=5, n.burnin=30, 
           model.file=model, jags.module=NULL)
mcmcfit = as.mcmc(fit)

ggsobject = ggs(mcmcfit)
ggs_density(ggsobject, "mu")
