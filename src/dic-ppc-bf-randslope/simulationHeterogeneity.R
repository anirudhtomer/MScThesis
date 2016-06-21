########## BAYESIAN PACKAGES ###########
install.packages('R2jags')
install.packages('ggmcmc')
install.packages("doParallel")
install.packages("MCMCpack")
install.packages("mvtnorm")
library(R2jags)
library(ggmcmc)
library(doParallel)
library(beepr)
library(MASS)
library(MCMCpack)
library(reshape2)
library(mvtnorm)
library(mixAK)

registerDoParallel(cores=8)

########## OTHER SOURCE CODE FILES ############
source("../common/fitModelRandSlopeConditional.R")
source("../common/generateDataRandSlope.R")
source("../common/createModelRandSlopeConditional.R")
source("../common/extractFuncRandSlopeConditional.R")
source("DIC_functions.R")
source("bfWishart.R")

########## Quick graphical analysis of the simulated mixture distribution #######
qplot(x=randIntercept, y=randSlope, data=data.frame(extractRandomComp(viaReg = T)), xlab="Random intercept", ylab="Random slope")

########## Running MCMC simulations ##########
numchains = 1
niter = 60000
nthin=80
nburnin=25000

ncomponents=3
if(ncomponents==1){
  fit = fitModel(niter, nthin, 0, jagsmodel = singleModel, nchains = numchains, ncomponents)  
  mcmcfit = as.mcmc(fit)
  colnames(mcmcfit[[1]])[match("Eta",colnames(mcmcfit[[1]]))] = "Eta[1]"
  colnames(mcmcfit[[1]])[match("precisionIntercept",colnames(mcmcfit[[1]]))] = "precisionIntercept[1]"
  colnames(mcmcfit[[1]])[match("precisionSlope",colnames(mcmcfit[[1]]))] = "precisionSlope[1]"
  colnames(mcmcfit[[1]])[match("rho",colnames(mcmcfit[[1]]))] = "rho[1]"
}else{
  fit = fitModel(niter, nthin, 0, jagsmodel = model, nchains = numchains, ncomponents)
  mcmcfit = as.mcmc(fit)
}

attributes(mcmcfit)$ncomponents = ncomponents
attributes(mcmcfit)$nsubjects = nsubjects
attributes(mcmcfit)$nrep = nrep
attributes(mcmcfit)$time = time
attributes(mcmcfit)$wishartPriorScale = wishartPriorScale

############## Calculating DIC #################
dic1=calculateDIC1(mcmcfit)
dic2=calculateDIC2(mcmcfit)
dic3=calculateDIC3(mcmcfit)
dic4=calculateDIC4(mcmcfit)
dic5=calculateDIC5(mcmcfit)
dic7=calculateDIC7(mcmcfit)

##### Calculating PPC ########
###### Run the code in ppc.R ###### 

#### Calculating marginal likelihood using Chib's approximation ####
#### Run bf.R one line and one function at a time ###

########## Graphical analysis of MCMC fit #########
ggs_density(ggsobject, "beta")
ggs_density(ggsobject, "errPrecision")
ggs_density(ggsobject, "randPrecision")
ggs_running(ggsobject, "randmu")
ggs_density(ggsobject, "Eta")

par(mfrow=c(2,2))
for(k in 1:attributes(mcmcfit)$ncomponents){
  for(i in 1:2){
    for(j in 1:2){
      paramName = paste("randSigma[",i,",",j,",",k,"]", sep="")
      parmData = mcmcfit[[1]][,paramName]
      plot(density(parmData), main=paramName)
      abline(v = median(parmData), col="red")
      abline(v = mean(parmData))
    }
  }
  readline()
}

ggs_compare_partial(ggsobject, "randmu", rug = T)

ggs_running(ggsobject, "beta")
ggs_running(ggsobject, "precisionIntercept")
ggs_running(ggsobject, "precisionSlope")
ggs_running(ggsobject, "rho")
ggs_running(ggsobject, "randmu")
ggs_running(ggsobject, "Eta")
ggs_running(ggsobject, "randSigma")

ggs_crosscorrelation(ggsobject)
ggs_pairs(ggsobject, "beta")
ggs_pairs(ggsobject, "Eta")
ggs_pairs(ggsobject, "randmu")
ggs_pairs(ggsobject, "Precision")

ggs_traceplot(ggsobject, "beta")
ggs_traceplot(ggsobject, "Precision")
ggs_traceplot(ggsobject, "randmu")
ggs_traceplot(ggsobject, "Eta")

########## Mode hunting: works only for 2 or 4 chains and when only random intercept is used in the model##########
mcmcfitdf = data.frame(y= numeric(0), x= numeric(0), 
              pairId = character(0), row = character(0), col=character(0))
rownum = c("1","1","2","2")
colnum = c("1","2","1","2")
for(i in 1:numchains){
  for(j in 1:ncomponents){
    for(k in 1:ncomponents){
      if(j!=k){
        y = as.numeric(mcmcfit[[i]][,paste("randmu[",k,"]", sep = "")])
        x = as.numeric(mcmcfit[[i]][,paste("randmu[",j,"]", sep = "")])
        pairId = rep(paste("mu(",sort(c(j,k)),")", sep=""), length(y))
        row=rep(rownum[i], length(y))
        col=rep(colnum[i], length(y))
        
        temp = data.frame(y,x,pairId,row, col)
        mcmcfitdf = rbind(mcmcfitdf, temp)
      }
    }
  }
}

xylims = c(100,200)
qplot(y=y, x=x, data=mcmcfitdf, color=pairId, facets = row~col, ylim = xylims, xlim=xylims) + 
  geom_abline(intercept = 0, slope = 1)
