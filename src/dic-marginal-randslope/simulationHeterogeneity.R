########## BAYESIAN PACKAGES ###########
install.packages('R2jags')
install.packages('mcmcse')
install.packages('ggmcmc')
library(R2jags)
library(ggmcmc)
library(doParallel)
library(beepr)

registerDoParallel(cores=8)

########## OTHER SOURCE CODE FILES ############
source("../common/fitModelRandSlope.R")
source("../common/generateDataRandSlope.R")
source("../common/createModelRandSlope.R")
source("../common/extractFuncRandSlope.R")
source("DIC_functions.R")

numchains = 1
niter = 3000
nthin=50
nburnin=200

ncomponents=3
fit = fitModel(niter, nthin, 0, jagsmodel = model, nchains = numchains, ncomponents)
mcmcfit = as.mcmc(fit)
attributes(mcmcfit)$ncomponents = ncomponents
attributes(mcmcfit)$nsubjects = nsubjects
attributes(mcmcfit)$nrep = nrep

beep(sound=8)

mcmcfit[[1]] = mcmcfit[[1]][((nburnin/nthin+1):(niter/nthin)),]
mcmcLen = nrow(mcmcfit[[1]])
attributes(mcmcfit[[1]])$mcpar = c(1, mcmcLen, 1)
ggsobject = ggs(mcmcfit)

calculateObsDIC1(mcmcfit, "obsDeviance")
calculateObsDIC2(mcmcfit, "obsDeviance")
calculateObsDIC3(mcmcfit, "obsDeviance","obsL")

dev.off()
########## Graphical analysis of the simulated mixture distribution #######
qplot(x=randIntercept, y=randSlope, data=data.frame(extractRandomComp(viaReg = T)))

########## Graphical analysis of MCMC fit #########
heidel.diag(mcmcfit)

ggs_density(ggsobject, "beta")
ggs_density(ggsobject, "Precision")
ggs_density(ggsobject, "randmu")
ggs_density(ggsobject, "Eta")

#Compare the whole chainw with the last part...similar to geweke
ggs_compare_partial(ggsobject, "randmu", rug = T)

ggs_running(ggsobject, "beta")
ggs_running(ggsobject, "Precision")
ggs_running(ggsobject, "randmu")
ggs_running(ggsobject, "Eta")

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

