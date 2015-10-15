install.packages("ggplot2")
library(ggplot2)

#mixture density
mat=rmultinom(1, 10000000, c(1/6,1/2,1/3))
first=rnorm(mat[1,1], -10, 3)
second=rnorm(mat[2,1], 0, 1)
third=rnorm(mat[3,1], 4, 2)
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

