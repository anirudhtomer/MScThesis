install.packages('R2jags')
install.packages('ggmcmc')
library(R2jags)
library(ggmcmc)

generateRandomSample = function(p1=-1,p2=-1,p3=-1,mu1,mu2,mu3){
  if(p1==-1 || p2==-1 || p3==-1){
    prob1=rbeta(n = 1, shape1 = 0.5,0.5)
    prob2=runif(1,0,1-prob1)
    mu1=rnorm(1,0,20)
    mu2=rnorm(1,0,20)
    mu3=rnorm(1,0,20)
  }else{
    prob1 = p1
    prob2 = p2
    prob3 = p3
  }
  
  mat=rmultinom(1, 800, c(prob1,prob2,1-(prob1+prob2)))
  mat=rmultinom(1, 800, c(prob1,prob2,1-(prob1+prob2)))
  first=rnorm(mat[1,1], mu1, 2)
  second=rnorm(mat[2,1], mu2, 2)
  third=rnorm(mat[3,1], mu3, 2)
  
  return (list("sample"=c(first, second, third), "eta"=c(prob1,prob2,1-prob1-prob2), "mu"=c(mu1,mu2,mu3)))
}

############ MIXTURE WITH INTERCEPT different etas ###############
randSamp=generateRandomSample()
sample = randSamp$sample
densityplot = ggplot()+ aes(sample) + geom_density()
densityplot + ylab(expression("p"[ Y ]*"(y)"))  + xlab("Y") + theme(axis.text=element_text(size=14),axis.title=element_text(size=18), plot.title=element_text(size=20))

intercept = 200
sample = sample + intercept

nobs = length(sample)
ncomponents=3
dirichParm = rep(1, ncomponents)

model=function(){
  for(i in 1:nobs){
    sample[i]~dnorm(mu[S[i]], precision)
    S[i]~dcat(eta[])
  }
  
  for(j in 1:ncomponents){
    mue[j]~dnorm(0, 0.0001)  
  }
  
  #for(m in 1:ncomponents){
  #  muk[m]<-mue[m]-(mue[1]*eta[1] + mue[2]*eta[2] + mue[3]*eta[3])
  #}
  mu<-sort(mue[])
  
  precision~dgamma(0.0005, 0.0005)
  intercept~dnorm(0, 0.0001)
  eta~ddirch(dirichParm[])
}

datanodes = c("sample","nobs","ncomponents", "dirichParm")
generateInitialValues=function(ncomp){
  list("precision"=c(1), 
       "eta"=rep(1/ncomponents,ncomponents),
       #"eta"=c(0.36,0.45,0.19),
       "mue"=rep(0,ncomponents),
       "intercept"=c(0))
}
initialValues = list(generateInitialValues(ncomponents),
                     generateInitialValues(ncomponents),
                     generateInitialValues(ncomponents),
                     generateInitialValues(ncomponents),
                     generateInitialValues(ncomponents))
params = c("precision","mu","eta","intercept")

unload.module("glm")
fit = jags(data=datanodes, inits=initialValues, params, 
           n.chains=length(initialValues), n.iter=7000,n.thin=50, n.burnin=1000, 
           model.file=model, jags.module=NULL)
mcmcfit = as.mcmc(fit)

ggsobject = ggs(mcmcfit)

ggs_density(ggsobject, "mu")
ggs_density(ggsobject, "intercept")
ggs_density(ggsobject, "eta")
ggs_running(ggsobject, "eta")
ggs_running(ggsobject, "mu")  
ggs_autocorrelation(ggsobject,"mu")
ggs_pairs(ggsobject)
