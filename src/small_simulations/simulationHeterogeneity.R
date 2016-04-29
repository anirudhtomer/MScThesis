########## BAYESIAN PACKAGES ###########
install.packages('R2jags')
install.packages('mcmcse')
install.packages('ggmcmc')
library(R2jags)
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

#make 4 groups
#weight is gender+by+quadtime+random intercept trimodal
#random intercept corresponds to diet---poor, ok, good
#time is 1 to 10
#males are 40 kg extra, 1997 is 30kg extra, quadtime coeff = 2, intercept is 50
#poor diet -5, ok diet 0, good diet +5

time = 1:10
weightgen=function(gender, by, diet){
  weight = rep(150, length(time));

  if(gender=="M"){
    weight = weight + 40
  }

  if(by=="97"){
    weight = weight + 30
  }

  weight = weight + time*10

  switch(diet,
         poor={
           weight = weight + rnorm(1, -30, 4)
         },
         ok={
           weight = weight + rnorm(1, 0, 4)
         },
         good={
           weight = weight + rnorm(1, 30, 4)
         })

  weight = weight + rnorm(length(time), 0, 3)
}

extractRandomComp = function(){
  temp = ds$weight
  temp = temp-150
  temp = temp-(as.numeric(ds$gender)-1)*40
  temp = temp-(as.numeric(ds$by)-1)*30
  temp = temp-ds$time*10
  temp
}

diets = c("poor", "ok", "good")
#60 males, 60 females, 30 each from 96 and 97
sublist = matrix(nrow = 0, ncol = 6)
for(gender in c("M", "F")){
  for(by in c("96", "97")){
    for(diet in diets){
      for(k in 1:15){
        id = rep(paste(gender, by, diet, k, sep=""), length(time))
        weight = weightgen(gender, by, diet)
        subject = cbind(id, rep(gender, length(time)),
                        rep(by, length(time)), time,
                        rep(diet, length(time)), weight)
        sublist=rbind(sublist, subject)
      }
    }
  }
}

ds = data.frame(sublist)
colnames(ds) = c("id", "gender", "by", "time", "diet", "weight")
ds$time =  as.numeric(as.character(ds$time))
ds$weight = as.numeric(as.character(ds$weight))
ds = ds[order(ds$id, ds$time),]
nobs = length(ds$id)
nfolks = nobs/length(time)
nrep = length(time)
ncomponents=length(diets)

extractRandomComp()

######### Bayesian hetergeneity model ########
betaPrecision = 0.00001
dirichParm = rep(1, ncomponents)

dsweight = ds$weight
dstime = ds$time
dsgender = as.numeric(ds$gender)-1
dsby = as.numeric(ds$by)-1

model=function(){
  for(i in 0:(nfolks-1)){
    for(j in 1:nrep){
      dsweight[i*nrep+j]~dnorm(mu[i*nrep+j],errPrecision)
      mu[i*nrep+j]<- betaGender*dsgender[i*nrep+j]+
        betaBy*dsby[i*nrep+j]+betaTime*dstime[i*nrep+j] + randomIntercept[i+1]
    }
    #randomIntercept[i+1]~dnorm(randmu[S[i+1]],randPrecision[S[i+1]])
    randomIntercept[i+1]~dnorm(randmu[S[i+1]],randPrecision)
    S[i+1]~dcat(eta[])
  }

  for(k in 1:ncomponents){
    #randPrecision[k]~dgamma(0.0005, 0.0005)
    randmu[k] ~ dnorm(0.0, betaPrecision)
  }
  randPrecision~dgamma(0.0005, 0.0005)

  #randmu[1:ncomponents]<-sort(randmue)

  errPrecision~dgamma(0.0005,0.0005)

  #Low precision is second parameter
  intercept~dnorm(0,betaPrecision)%_%I(149.8,151.2)
  betaGender~dnorm(0,betaPrecision)
  betaBy~dnorm(0,betaPrecision)
  betaTime~dnorm(0,betaPrecision)

  eta~ddirch(dirichParm[])
}

initS=function(nfolks, ncomponents){
  retval = numeric(nfolks)
  randomorder=sample(1:nfolks, nfolks)
  for(i in 1:ncomponents){
    retval[randomorder[((i-1)*(nfolks/ncomponents)+1):(i*(nfolks/ncomponents))]] = i;
  }
  return (retval)
}

datanodes = c("dsweight","dsby","dsgender","dstime","nfolks", "nrep", "betaPrecision", "ncomponents", "dirichParm")
initialValues = list(list("betaGender"=c(0),"betaBy"=c(0),
                          "betaTime"=c(0),"errPrecision"=c(1),"intercept"=c(0),
                          "randPrecision"=rep(1, 1), "eta"=rep(1/ncomponents, ncomponents),
                          "randmu"=rep(0, ncomponents), "S"=initS(nfolks, ncomponents)))
stochasticNodes = c("intercept","betaGender","betaBy","betaTime","errPrecision", "randPrecision", "randmu", "eta")

unload.module("glm")
fit = jags(data=datanodes, inits=initialValues, stochasticNodes,
           n.chains=1, n.iter=50000,n.thin=50, n.burnin=500,
           model.file=model, jags.module=NULL)
mcmcfit = as.mcmc(fit)

#Resulting values
for(i in 1:4){
  print(paste(stochasticNodes[i], mean(mcmcfit[,stochasticNodes[i]])))
  plotdensity(mcmcfit[,stochasticNodes[i]], stochasticNodes[i])
}
print(paste("errPrecision", sqrt(1/mean(mcmcfit[,"errPrecision"]))))
for(i in 6:8){
  for(j in 1:ncomponents){
    nodename = paste(stochasticNodes[i],"[",j,"]", sep = "")
    if(grepl("Preci+", i)){
      print(paste(nodename, sqrt(1/mean(mcmcfit[,nodename]))))
    }else{
      print(paste(nodename, mean(mcmcfit[,nodename])))
      plotdensity(mcmcfit[,nodename], nodename)
    }
  }
}
