install.packages('R2OpenBUGS')
install.packages("mcmcplots")
install.packages("superdiag")
install.packages('mcmcse')
install.packages('ggmcmc')
library(R2OpenBUGS)
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

############ MIXTURE ###############
mat=rmultinom(1, 1000, c(1/2,1/2))
first=rnorm(mat[1,1], 10, 8)
second=rnorm(mat[2,1], -10, 8)
sample = c(first, second)

densityplot = ggplot()+ aes(sample) + geom_density()
densityplot + ylab(expression("p"[ Y ]*"(y)"))  + xlab("Y")

nobs = length(sample)
ncomponents=2
dirichParm = rep(1, ncomponents)

model=function(){
  for(i in 1:nobs){
    sample[i]~dnorm(mu[S[i]], precision[S[i]])
    S[i]~dcat(eta[])
  }
  
  for(j in 1:ncomponents){
    precision[j]~dgamma(0.005, 0.005)
  }
  
  mu[1]~dnorm(0, 0.0001)
  mu[2]~dnorm(0, 0.0001)
  
  eta[1:ncomponents]~ddirich(dirichParm[])
}

datanodes = list("sample"=sample,"nobs"=nobs,"ncomponents"=ncomponents, 
                 "dirichParm"=dirichParm)
initialValues = function(){
  list("eta"=rep(1/ncomponents, ncomponents),
       "mu"=rep(0,ncomponents), "precision"=rep(1,ncomponents))
}
params = c("precision","mu","eta")

fit = bugs(data=datanodes, inits=initialValues, parameters.to.save = params, 
           n.chains=1, n.iter=10000,n.thin=5, n.burnin=30, 
           model.file=model)
mcmcfit = as.mcmc(fit)

plotdensity(mcmcfit[,"precision[1]"], "precision[1]")
plotdensity(mcmcfit[,"precision[2]"], "precision[1]")
plotdensity(mcmcfit[,"mu[1]"], "mu[1]")
plotdensity(mcmcfit[,"mu[2]"], "mu[2]")
plotdensity(mcmcfit[,"eta[1]"], "eta[1]")
plotdensity(mcmcfit[,"eta[2]"], "eta[2]")

############ MIXTURE WITH INTERCEPT ###############
mat=rmultinom(1, 1000, c(1/2,1/2))
first=rnorm(mat[1,1], 10, 8)
second=rnorm(mat[2,1], -10, 8)
sample = c(first, second)

sample = sample + 40

densityplot = ggplot()+ aes(sample) + geom_density()
densityplot + ylab(expression("p"[ Y ]*"(y)"))  + xlab("Y")

nobs = length(sample)
ncomponents=2
dirichParm = rep(1, ncomponents)

model=function(){
  for(i in 1:nobs){
    sample[i]~dnorm(temp[i], precision[S[i]])
    temp[i]<-intercept + mu[S[i]]
    S[i]~dcat(eta[])
  }
  
  for(j in 1:ncomponents){
    precision[j]~dgamma(0.001, 0.001)
  }
  
  intercept~dnorm(0, 0.0001)%_%I(39.5,40.5)
  mu[1]~dnorm(0,0.0001)
  mu[2]~dnorm(0,0.0001)
  
  eta[1:ncomponents]~ddirich(dirichParm[])
}

datanodes = list("sample"=sample,"nobs"=nobs,"ncomponents"=ncomponents, 
                 "dirichParm"=dirichParm)
initialValues = function(){
  list("eta"=rep(1/ncomponents, ncomponents),
       "mu"=rep(0,ncomponents),"intercept"=0, "precision"=rep(1,ncomponents))
}
params = c("precision","mu","eta","intercept")

fit = bugs(data=datanodes, inits=initialValues, parameters.to.save = params, 
           n.chains=1, n.iter=10000,n.thin=5, n.burnin=30, 
           model.file=model)
mcmcfit = as.mcmc(fit)

plotdensity(mcmcfit[,"precision[1]"], "precision[1]")
plotdensity(mcmcfit[,"precision[2]"], "precision[1]")
plotdensity(mcmcfit[,"mu[1]"], "mu[1]")
plotdensity(mcmcfit[,"mu[2]"], "mu[2]")
plotdensity(mcmcfit[,"intercept"], "intercept")
plotdensity(mcmcfit[,"eta[1]"], "eta[1]")
plotdensity(mcmcfit[,"eta[2]"], "eta[2]")

############ MIXTURE WITH RANDOM ERROR EVERY 10 OBS###############
#######data extracted from another program sorry ##########
sample = extractRandomComp()

densityplot = ggplot()+ aes(sample) + geom_density()
densityplot + ylab(expression("p"[ Y ]*"(y)"))  + xlab("Y")

nobs = length(sample)
nrep = 10
nfolks = nobs/nrep
ncomponents=2
dirichParm = rep(1, ncomponents)

model=function(){
  for(i in 0:(nfolks-1)){
    for(j in 1:nrep){
      sample[i*nrep+j]~dnorm(randomIntercept[i+1],errPrecision)
    }
    randomIntercept[i+1]~dnorm(randmu[S[i+1]],randPrecision[S[i+1]])
    S[i+1]~dcat(eta[])
  }
  
  errPrecision~dgamma(0.005,0.005)
  
  for(j in 1:ncomponents){
    randPrecision[j]~dgamma(0.005, 0.005)
  }
  
  randmu[1]~dnorm(0, 0.0001)
  randmu[2]~dnorm(0, 0.0001)
  
  eta[1:ncomponents]~ddirich(dirichParm[])
}

datanodes = list("sample"=sample,"nrep"=nrep,"ncomponents"=ncomponents, 
                 "dirichParm"=dirichParm, "nfolks"=nfolks)
initialValues = function(){
  list("eta"=rep(1/ncomponents, ncomponents),"errPrecision"=1,
       "randmu"=rep(0,ncomponents), "randPrecision"=rep(1,ncomponents))
}
params = c("randPrecision","errPrecision","randmu","eta")

fit = bugs(data=datanodes, inits=initialValues, parameters.to.save = params, 
           n.chains=1, n.iter=10000,n.thin=5, n.burnin=30, 
           model.file=model)
mcmcfit = as.mcmc(fit)

plotdensity(mcmcfit[,"precision[1]"], "precision[1]")
plotdensity(mcmcfit[,"precision[2]"], "precision[1]")
plotdensity(mcmcfit[,"mu[1]"], "mu[1]")
plotdensity(mcmcfit[,"mu[2]"], "mu[2]")
plotdensity(mcmcfit[,"eta[1]"], "eta[1]")
plotdensity(mcmcfit[,"eta[2]"], "eta[2]")
