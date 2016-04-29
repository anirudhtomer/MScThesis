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

numchains = 1
niter = 80000
nthin=50
nburnin=30000

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

#############################

for (ncomponents in 2:8){
fit = fitModel(niter, nthin, 0, jagsmodel = model, nchains = numchains, ncomponents)
mcmcfit = as.mcmc(fit)
attributes(mcmcfit)$ncomponents = ncomponents
attributes(mcmcfit)$nsubjects = nsubjects
attributes(mcmcfit)$nrep = nrep
attributes(mcmcfit)$time = time
attributes(mcmcfit)$wishartPriorScale = wishartPriorScale

if(ncomponents==2){
  mcmcfit_5fused2comp=mcmcfit
}

if(ncomponents==3){
  mcmcfit_5fused3comp = mcmcfit 
}

if(ncomponents==4){
  mcmcfit_5fused4comp= mcmcfit 
}

if(ncomponents==5){
  mcmcfit_5fused5comp=mcmcfit
}

if(ncomponents==6){
  mcmcfit_5fused6comp=mcmcfit
}

if(ncomponents==7){
  mcmcfit_5fused7comp=mcmcfit
}
if(ncomponents==8){
  mcmcfit_5fused8comp=mcmcfit
}
  
save.image("D:/Dropbox/MSc Stats/Thesis/MScThesis/latex/code/dic-ppc-bf-randslope/wishart.RData")
}

pd=numeric()
dic = numeric()

mcmcfit[[1]] = mcmcfit[[1]][((nburnin/nthin+1):(niter/nthin)),]
mcmcLen = nrow(mcmcfit[[1]])
attributes(mcmcfit[[1]])$mcpar = c(1, mcmcLen, 1)
nrow(mcmcfit[[1]])

dic1=calculateDIC1(mcmcfit)
dic2=calculateDIC2(mcmcfit)
dic3=calculateDIC3(mcmcfit)
dic4=calculateDIC4(mcmcfit)
dic5=calculateDIC5(mcmcfit)
dic7=calculateDIC7(mcmcfit)

sapply(list(dic1,dic2,dic3,dic4,dic5,dic7), FUN = function(x){pd[length(pd)+1]<<-x$pd; dic[length(dic)+1]<<-x$dic}) 
write.csv(x=data.frame(pd, dic), file = "C:/Users/Anirudh/Desktop/temp.csv", row.names = F)


##############

beep(sound=8)

mcmcfit[[1]] = mcmcfit[[1]][((nburnin/nthin+1):(niter/nthin)),]
mcmcLen = nrow(mcmcfit[[1]])
attributes(mcmcfit[[1]])$mcpar = c(1, mcmcLen, 1)
ggsobject = ggs(mcmcfit)

attributes(mcmcfit)
#save.image("F:/docs/Dropbox/MSc Stats/Thesis/MScThesis/latex/code/dic-ppc-randslope/.RData")

########## Graphical analysis of the simulated mixture distribution #######
qplot(x=randIntercept, y=randSlope, data=data.frame(extractRandomComp(viaReg = T)))

########## Graphical analysis of MCMC fit #########
heidel.diag(mcmcfit)

ggs_density(ggsobject, "beta")
ggs_density(ggsobject, "errPrecision")
ggs_density(ggsobject, "randSigma")
ggs_density(ggsobject, "randmu")
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

#Compare the whole chain with the last part...similar to geweke
ggs_compare_partial(ggsobject, "randmu", rug = T)

ggs_running(ggsobject, "beta")
ggs_running(ggsobject, "precisionIntercept")
ggs_running(ggsobject, "precisionSlope")
ggs_running(ggsobject, "rho")
ggs_running(ggsobject, "randmu")
ggs_running(ggsobject, "Eta")
ggs_running(ggsobject, "randSigma")

install.packages("igraph")
library(igraph)
plot(running.mean(mcmcfit[[1]][,"randSigma[1,1,1]"], binwidth = 2))

ggs_autocorrelation(ggsobject, "beta")
ggs_autocorrelation(ggsobject, "Precision")
ggs_autocorrelation(ggsobject, "randmu")
ggs_autocorrelation(ggsobject, "Eta")

#dont attempt ggs_pairs on entire object
ggs_crosscorrelation(ggsobject)
ggs_pairs(ggsobject, "beta")
ggs_pairs(ggsobject, "Eta")
ggs_pairs(ggsobject, "randmu")
ggs_pairs(ggsobject, "Precision")

ggs_traceplot(ggsobject, "beta")
ggs_traceplot(ggsobject, "Precision")
ggs_traceplot(ggsobject, "randmu")
ggs_traceplot(ggsobject, "Eta")

########## Mode hunting, works only for 2 or 4 chains ##########
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
