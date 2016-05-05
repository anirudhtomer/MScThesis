obsCovMatrix =  var(splaash1)
priorDf = 3
posteriorScale = solve(solve(diag(2)) + obsCovMatrix*59)
posteriordf = priorDf + 60

var1 = numeric()
var2 = numeric()
cov = numeric()

for(i in 1:1000){
  posteriorPrecisionSample = rwish(wishartPostDf[1], wishartPosteriorScale[[1]])
  posteriorVarianceSample  = (posteriorPrecisionSample)
  var1[i] = posteriorVarianceSample[1,1]
  var2[i] = posteriorVarianceSample[2,2]
  cov[i] = posteriorVarianceSample[2,1]
} 

par(mfrow=c(1,3))
plot(density(var1))
plot(density(var2))
plot(density(cov))

meanSplaash1 = c(mean(splaash1[,1]),mean(splaash1[,2]))

corrrr=c()
for(i in 1:500){
  temp=solve(rWishart(n=1,df=2, diag(2))[,,1])
  corrrr[i] = temp[1,1]
}
plot(density(corrrr))