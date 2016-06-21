########## BAYESIAN PACKAGES ###########
install.packages('R2jags')
install.packages('ggmcmc')
install.packages("doParallel")
library(R2jags)
library(ggmcmc)
library(doParallel)
library(beepr)
library(MASS)
library(reshape2)
library(corrplot)

registerDoParallel(cores=8)

########## OTHER SOURCE CODE FILES ############
source("datacleaning.R")

source("fitModel.R")
source("createModel.R")
source("../common/extractFuncRandSlopeConditional.R")
source("DIC_functions.R")
source("bf.R")

########## Quick graphical analysis of the random component of the model #######
qplot(x=randIntercept, y=randSlope, data=data.frame(extractRandomComp()), xlim = c(-25,25),ylim = c(-25,25))
qplot(y=Hb, x=TSPD,data=ds[1:400,], geom = "line", group=Id)

numchains= 1
niter = 1000000
nthin=200
nburnin=700000

ncomponents=3
if(ncomponents==1){
  fit = fitModel(niter, nthin, 0, jagsmodel = singleModel, nchains = numchains, ncomponents, 1)  
  mcmcfit = as.mcmc(fit)
  colnames(mcmcfit[[1]])[match("Eta",colnames(mcmcfit[[1]]))] = "Eta[1]"
  colnames(mcmcfit[[1]])[match("precision1",colnames(mcmcfit[[1]]))] = "precision1[1]"
  colnames(mcmcfit[[1]])[match("precision2",colnames(mcmcfit[[1]]))] = "precision2[1]"
  colnames(mcmcfit[[1]])[match("rho",colnames(mcmcfit[[1]]))] = "rho[1]"
}else{
  fit = fitModel(niter, nthin, 0, jagsmodel = model, nchains = numchains, ncomponents, 2.3)
  mcmcfit = as.mcmc(fit)
}

attributes(mcmcfit)$ncomponents = ncomponents
attributes(mcmcfit)$nsubjects = nsubjects
getRandSigma(mcmcfit)


calculateDIC1(mcmcfit)
calculateDIC2(mcmcfit)
calculateDIC3(mcmcfit)
calculateDIC4(mcmcfit)
calculateDIC5(mcmcfit)
calculateDIC7(mcmcfit)

####### for ppc run the code in ppc.R########