install.packages("ggplot2")
library(ggplot2)

#mixture density
mat=rmultinom(1, 1000, c(1/6,1/2,1/3))
first=rnorm(mat[1,1], -10, 3)
second=rnorm(mat[2,1], 0, 1)
third=rnorm(mat[3,1], 10, 3)
sample = c(first, second, third)

densityplot = ggplot()+ aes(sample) + geom_density()
densityplot + ylab(expression("p"[ Y ]*"(y)"))  + xlab("Y") + theme(axis.text=element_text(size=14),axis.title=element_text(size=18), plot.title=element_text(size=20))

#toy problem-beta prior
sample = rbeta(n=100000, shape1 = 9, shape2 = 3)
densityplot = ggplot()+ aes(sample) + geom_density()
densityplot + ylab(expression("p("*pi*")"))  + xlab(expression(pi)) + theme(axis.text=element_text(size=14),axis.title=element_text(size=18), plot.title=element_text(size=20))

#toy problem-beta posterior
sample = rbeta(n=100000, shape1 = 9 + 6, shape2 = 3 + 10 - 6)
densityplot = ggplot()+ aes(sample) + geom_density()
densityplot + ylab(expression("p("*pi*")"))  + xlab(expression(pi)) + theme(axis.text=element_text(size=14),axis.title=element_text(size=18), plot.title=element_text(size=20))

#example for motivation of mixture random effects
group1 = 50
group2 = 90

subject_count = 30
variance1=10
variance2=10

response1 = c(rnorm(subject_count, group1, variance1), rnorm(subject_count, group2, variance2))
response2 = c(rnorm(subject_count, group1, variance1) * 2, rnorm(subject_count, group2, variance2)*2)
response3 = c(rnorm(subject_count, group1, variance1) * 2.5, rnorm(subject_count, group2, variance2)*2.5)
response4 = c(rnorm(subject_count, group1, variance1) * 3, rnorm(subject_count, group2, variance2)*3)
response5 = c(rnorm(subject_count, group1, variance1) * 3.5, rnorm(subject_count, group2, variance2)*3.5)
response6 = c(rnorm(subject_count, group1, variance1) * 4, rnorm(subject_count, group2, variance2)*4)

resp = data.frame(response1, response2, response3, response4, response5, response6)
plot(NULL,xlab="Time", ylab = "Response (Y)", xlim=c(0, 5), ylim=c(0, 400),type = 'n', yaxs="i", xaxs="i")

for(i in 1:(subject_count*2)){
  if(i>subject_count){
    color="blue"
  }else{
    color="darkgreen"
  }
  lines(x = c(0:5), y=c(resp$response1[i], resp$response2[i], resp$response3[i], resp$response4[i], resp$response5[i], resp$response6[i]), col=color)
}
legend("topleft", lty=1, col=c("blue", "darkgreen"),  legend=c("High risk group", "Low risk group"))


#######  Constrained jags trial ########

model=function(){
  for(i in 1:10){
    y1[i]~dnorm(mu1,1)
    y2[i]~dnorm(mu2,1)
  }
  mu1~dnorm(0,1)
  z~dinterval(mu2*mu2, cons-mu1*mu1)
  mu2~dnorm(0,1)
  
  #z~dbern(constraint)
  #constraint <- step(mu1*mu1 + mu2*mu2 - cons)
}

y1 = rnorm(10, 0, 1)
y2 = rnorm(10, 0, 1)
y1 = c(-0.35, -0.46,  1.05, -0.20, -0.70,  0.19,  1.07,  1.10,  0.54, -0.01)
y2 = c(0.09, 0.25,  0.04, -1.03,  0.53,  1.56,  0.03, -0.13,  0.31,  0.54)
z=1

cons = 0.2
 
datanodes = list("y1"=y1,"y2"=y2,"z"=z,"cons"=cons)
initialValues = list(list("mu1"=c(0.5), "mu2"=c(0.5)))
stochasticNodes = c("mu1", "mu2", "z")
 
unload.module("glm")
fit = jags(data=datanodes, inits=initialValues, 
           parameters.to.save = stochasticNodes, n.chains=1, 
           n.iter=20000, n.thin=5, n.burnin=50, model.file=model, 
           jags.module=NULL)
mcmcfit = as.mcmc(fit)
ggs_density(ggs(mcmcfit))
mcmcdf = data.frame(mcmcfit[[1]])
mcmcdf$cons = (mcmcdf$mu1^2 + mcmcdf$mu2^2)>cons
qplot(y=mu1,x=mu2, data=mcmcdf, color=cons)
