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

plotdensity=function(datasample, title="",ymax=-1){
  temp=density(datasample)
  if(ymax==-1){
    plot(temp, main=title, type="n")  
  }else{
    plot(temp, main=title, type="n", ylim = c(0,ymax))
  }
  
  modeIndex = order(temp$y, decreasing = T)[1]
  polygon(temp, col="lightgray", border="gray")
  segments(x0 = temp$x[modeIndex],y0=0, x1=temp$x[modeIndex],y1=temp$y[modeIndex],col = "red")
}

#mixture density
mat=rmultinom(1, 1000, c(1/6,1/2,1/3))
first=rnorm(mat[1,1], -10, 3)
second=rnorm(mat[2,1], 0, 1)
third=rnorm(mat[3,1], 10, 3)
sample = c(first, second, third)

nobs = length(sample)
ncomponents=3
dirichParm = rep(1, ncomponents)

model=function(){
  for(i in 1:nobs){
    sample[i]~dnorm(mu[S[i]], precision[S[i]])
    S[i]~dcat(eta[])
  }
  
  for(j in 1:ncomponents){
    mu[j]~dnorm(0, 0.00001)  
    precision[j]~dgamma(0.0005, 0.0005)
  }
  
  eta~ddirch(dirichParm[])
}

datanodes = c("sample","nobs","ncomponents", "dirichParm")
initialValues = list(list("precision"=rep(1, ncomponents), 
                          "eta"=rep(1/ncomponents, ncomponents),
                          "mu"=rep(0, ncomponents)))
params = c("precision","mu","eta","S")

unload.module("glm")
fit = jags(data=datanodes, inits=initialValues, params, 
           n.chains=1, n.iter=30000,n.thin=5, n.burnin=30, 
           model.file=model, jags.module=NULL)
mcmcfit = as.mcmc(fit)

plotdensity(mcmcfit[,"precision[1]"], "precision[1]")
plotdensity(mcmcfit[,"precision[2]"], "precision[1]")
plotdensity(mcmcfit[,"precision[3]"], "precision[1]")
plotdensity(mcmcfit[,"mu[1]"], "mu[1]")
plotdensity(mcmcfit[,"mu[2]"], "mu[2]")
plotdensity(mcmcfit[,"mu[3]"], "mu[3]")
plotdensity(mcmcfit[,"eta[1]"], "eta[1]")
plotdensity(mcmcfit[,"eta[2]"], "eta[2]")
plotdensity(mcmcfit[,"eta[3]"], "eta[3]")

############ MIXTURE WITH INTERCEPT ###############

mat=rmultinom(1, 1000, c(1/6,1/2,1/3))
first=rnorm(mat[1,1], -16, 3)
second=rnorm(mat[2,1], 0, 3)
third=rnorm(mat[3,1], 16, 3)
sample = c(first, second, third)

densityplot = ggplot()+ aes(sample) + geom_density()
densityplot + ylab(expression("p"[ Y ]*"(y)"))  + xlab("Y") + theme(axis.text=element_text(size=14),axis.title=element_text(size=18), plot.title=element_text(size=20))

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
  mu[1:ncomponents]<-sort(mue)
  
  precision~dgamma(0.0005, 0.0005)
  intercept~dnorm(0, 0.0001)%_%I(18.9,20.1)
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

plotdensity(mcmcfit[,"precision"], "precision")
plotdensity(mcmcfit[,"mu[1]"], "mu[1]")
plotdensity(mcmcfit[,"mu[2]"], "mu[2]")
plotdensity(mcmcfit[,"mu[3]"], "mu[3]")
plotdensity(mcmcfit[,"eta[1]"], "eta[1]")
plotdensity(mcmcfit[,"eta[2]"], "eta[2]")
plotdensity(mcmcfit[,"eta[3]"], "eta[3]")
plotdensity(mcmcfit[,"intercept"], "intercept")


