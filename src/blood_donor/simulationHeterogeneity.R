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

beep(sound=8)

mcmcfit[[1]] = mcmcfit[[1]][((nburnin/nthin+1):(niter/nthin)),]
mcmcLen = nrow(mcmcfit[[1]])
attributes(mcmcfit[[1]])$mcpar = c(1, mcmcLen, 1)
ggsobject = ggs(mcmcfit)

attributes(mcmcfit)

ncomponents=1
fit = fitModel(niter, nthin, 0, jagsmodel = singleModel, nchains = numchains, ncomponents, 1)  
mcmcfit_1comp_morecov = as.mcmc(fit)
attributes(mcmcfit_1comp_morecov)$ncomponents = 1
attributes(mcmcfit_1comp_morecov)$nsubjects = nsubjects
attributes(mcmcfit_1comp_morecov)$INTERCEPT_SCALE = INTERCEPT_SCALE
colnames(mcmcfit_1comp_morecov[[1]])[match("Eta",colnames(mcmcfit_1comp_morecov[[1]]))] = "Eta[1]"
colnames(mcmcfit_1comp_morecov[[1]])[match("precision1",colnames(mcmcfit_1comp_morecov[[1]]))] = "precision1[1]"
colnames(mcmcfit_1comp_morecov[[1]])[match("precision2",colnames(mcmcfit_1comp_morecov[[1]]))] = "precision2[1]"
colnames(mcmcfit_1comp_morecov[[1]])[match("rho",colnames(mcmcfit_1comp_morecov[[1]]))] = "rho[1]"
save.image("D:/Dropbox/MSc Stats/Thesis/MScThesis/src/blood_donor/.RData")

ncomponents=2
fit = fitModel(niter, nthin, 0, jagsmodel = model, nchains = numchains, ncomponents, 2.7)
mcmcfit_2comp_dir2dot7 = as.mcmc(fit)
attributes(mcmcfit_2comp_dir2dot7)$ncomponents = 2
attributes(mcmcfit_2comp_dir2dot7)$nsubjects = nsubjects
attributes(mcmcfit_2comp_dir2dot7)$INTERCEPT_SCALE = INTERCEPT_SCALE
save.image("D:/Dropbox/MSc Stats/Thesis/MScThesis/src/blood_donor/.RData")

ncomponents=3
fit = fitModel(niter, nthin, 0, jagsmodel = model, nchains = numchains, ncomponents, 0.5)
mcmcfit_3comp = as.mcmc(fit)
attributes(mcmcfit_3comp)$ncomponents = 3
attributes(mcmcfit_3comp)$nsubjects = nsubjects
attributes(mcmcfit_3comp)$INTERCEPT_SCALE = INTERCEPT_SCALE
save.image("D:/Dropbox/MSc Stats/Thesis/MScThesis/src/blood_donor/.RData")

ncomponents=4
fit = fitModel(niter, nthin, 0, jagsmodel = model, nchains = numchains, ncomponents, 0.5)
mcmcfit_4comp = as.mcmc(fit)
attributes(mcmcfit_4comp)$ncomponents = 4
attributes(mcmcfit_4comp)$nsubjects = nsubjects
attributes(mcmcfit_4comp)$INTERCEPT_SCALE = INTERCEPT_SCALE
save.image("D:/Dropbox/MSc Stats/Thesis/MScThesis/src/blood_donor/.RData")

ncomponents=5
fit = fitModel(niter, nthin, 0, jagsmodel = model, nchains = numchains, ncomponents, 3)
mcmcfit_5comp_dirich1 = as.mcmc(fit)
attributes(mcmcfit_5comp_dirich1)$ncomponents = 5
attributes(mcmcfit_5comp_dirich1)$nsubjects = nsubjects
attributes(mcmcfit_5comp_dirich1)$INTERCEPT_SCALE = INTERCEPT_SCALE
save.image("D:/Dropbox/MSc Stats/Thesis/MScThesis/src/blood_donor/.RData")

ncomponents=6
fit = fitModel(niter, nthin, 0, jagsmodel = model, nchains = numchains, ncomponents, 2)
mcmcfit_6comp = as.mcmc(fit)
attributes(mcmcfit_6comp)$ncomponents = 6
attributes(mcmcfit_6comp)$nsubjects = nsubjects
attributes(mcmcfit_6comp)$INTERCEPT_SCALE = INTERCEPT_SCALE
save.image("D:/Dropbox/MSc Stats/Thesis/MScThesis/src/blood_donor/.RData")

mcmcfit = mcmcfit_4comp
mcmcfit[[1]] = mcmcfit[[1]][((nburnin/nthin+1):(niter/nthin)),]
mcmcLen = nrow(mcmcfit[[1]])
attributes(mcmcfit[[1]])$mcpar = c(1, mcmcLen, 1)
dic4comp = list(calculateDIC1(mcmcfit), calculateDIC2(mcmcfit), calculateDIC3(mcmcfit), calculateDIC4(mcmcfit), calculateDIC5(mcmcfit), calculateDIC7(mcmcfit))
save.image("D:/Dropbox/MSc Stats/Thesis/MScThesis/src/blood_donor/.RData")

mcmcfit = mcmcfit_5comp_dirich1
mcmcfit[[1]] = mcmcfit[[1]][((nburnin/nthin+1):(niter/nthin)),]
mcmcLen = nrow(mcmcfit[[1]])
attributes(mcmcfit[[1]])$mcpar = c(1, mcmcLen, 1)
dic5comp = list(calculateDIC1(mcmcfit), calculateDIC2(mcmcfit), calculateDIC3(mcmcfit), calculateDIC4(mcmcfit), calculateDIC5(mcmcfit), calculateDIC7(mcmcfit))
save.image("D:/Dropbox/MSc Stats/Thesis/MScThesis/src/blood_donor/.RData")

mcmcfit = mcmcfit_6comp
mcmcfit[[1]] = mcmcfit[[1]][((nburnin/nthin+1):(niter/nthin)),]
mcmcLen = nrow(mcmcfit[[1]])
attributes(mcmcfit[[1]])$mcpar = c(1, mcmcLen, 1)
dic6comp = list(calculateDIC1(mcmcfit), calculateDIC2(mcmcfit), calculateDIC3(mcmcfit), calculateDIC4(mcmcfit), calculateDIC5(mcmcfit), calculateDIC7(mcmcfit))
save.image("D:/Dropbox/MSc Stats/Thesis/MScThesis/src/blood_donor/.RData")

calculateDIC1(mcmcfit)
calculateDIC2(mcmcfit)
calculateDIC3(mcmcfit)
calculateDIC4(mcmcfit)
calculateDIC5(mcmcfit)
calculateDIC7(mcmcfit)

########## Graphical analysis of the simulated mixture distribution #######
qplot(x=randIntercept, y=randSlope, data=data.frame(extractRandomComp()), xlim = c(-25,25),ylim = c(-25,25))
qplot(y=Hb, x=TSPD,data=ds[1:400,], geom = "line", group=Id)

########## Graphical analysis of MCMC fit #########
heidel.diag(mcmcfit)

ggsobject = ggs(mcmcfit_3comp)
ggs_density(ggsobject, "beta")
ggs_density(ggsobject, "errPrecision")
ggs_density(ggsobject, "randSigma")
ggs_density(ggsobject, "randmu")
ggs_density(ggsobject, "Eta")

par(mfrow=c(2,2))
for(k in 1:ncomponents){
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

install.packages("latex2exp")
library(latex2exp)

qplot(mcmcfit[[1]][, "randmu[1,1]"], geom="density", xlab=TeX("b_{Intercept-1}^C | y"), ylab="Density")
qplot(mcmcfit[[1]][, "randmu[1,2]"], geom="density", xlab=TeX("b_{Slope-1}^C | y"), ylab="Density")
qplot(mcmcfit[[1]][, "randmu[2,1]"], geom="density", xlab=TeX("b_{Intercept-2}^C | y"), ylab="Density")
qplot(mcmcfit[[1]][, "randmu[2,2]"], geom="density", xlab=TeX("b_{Slope-2}^C | y"), ylab="Density")

qplot(mcmcfit[[1]][, "Eta[1]"], geom="density", xlab=TeX("$\\eta_1 | y$"), ylab="Density")
qplot(mcmcfit[[1]][, "Eta[2]"], geom="density", xlab=TeX("$\\eta_2 | y$"), ylab="Density")

qplot(mcmcfit[[1]][, "randSigma[1,1,1]"], geom="density", xlab=TeX("$\\sigma^2_{Intercept-1}$"), ylab="Density")
qplot(mcmcfit[[1]][, "randSigma[1,2,1]"], geom="density", xlab=TeX("$\\sigma_{Intercept-1, Slope-1}$"), ylab="Density")
qplot(mcmcfit[[1]][, "randSigma[2,2,1]"], geom="density", xlab=TeX("$\\sigma^2_{Slope-1}$"), ylab="Density")

qplot(mcmcfit[[1]][, "randSigma[1,1,2]"], geom="density", xlab=TeX("$\\sigma^2_{Intercept-2}$"), ylab="Density")
qplot(mcmcfit[[1]][, "randSigma[1,2,2]"], geom="density", xlab=TeX("$\\sigma_{Intercept-2, Slope-2}$"), ylab="Density")
qplot(mcmcfit[[1]][, "randSigma[2,2,2]"], geom="density",  xlab=TeX("$\\sigma^2_{Slope-2}$"), ylab="Density")

#Compare the whole chain with the last part...similar to geweke
ggs_compare_partial(ggsobject, "randmu")

ggs_running(ggsobject, "beta")
ggs_running(ggsobject, "Precision")
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

temp=function(ddd){
  m1 = round(mean(ddd),3); 
  med1 = round(median(ddd),3); 
  med3 = c(round(HPDinterval(mcmc(ddd)),3))
  
  print(paste(" &",m1,"&", med1, "&", med3[1], ",", med3[2], sep=" "))
}
