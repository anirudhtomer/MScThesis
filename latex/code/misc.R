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

#toy problem 
sample = rbeta(n=100000, shape1 = 9, shape2 = 3)
densityplot = ggplot()+ aes(sample) + geom_density()
densityplot + ylab(expression("p("*pi*")"))  + xlab(expression(pi)) + theme(axis.text=element_text(size=14),axis.title=element_text(size=18), plot.title=element_text(size=20))


