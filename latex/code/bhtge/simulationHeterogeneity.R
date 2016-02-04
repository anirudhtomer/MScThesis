########## BAYESIAN PACKAGES ###########
install.packages('R2jags')
install.packages('mcmcse')
install.packages('ggmcmc')
library(R2jags)
library(mcmcse)
library(ggmcmc)

########## OTHER SOURCE CODE FILES ############
source("fitModel.R")
source("generateData.R")
source("createModel.R")
source("plotMCMC.R")

numchains = 5
fit = fitModel(niter = 20000, jagsmodel = model, nchains = numchains)
mcmcfit = as.mcmc(fit)
ggsobject = ggs(mcmcfit)

dev.off()

heidel.diag(mcmcfit)

ggs_density(ggsobject, "beta", rug=T)
ggs_density(ggsobject, "Precision", rug=T)
ggs_density(ggsobject, "randmu", rug=T)
ggs_density(ggsobject, "Eta",rug=T)

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
