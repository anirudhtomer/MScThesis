install.packages("ggplot2")
library(ggplot2)

mat=rmultinom(1, 10000000, c(1/6,1/2,1/3))
first=rnorm(mat[1,1], -10, 3)
second=rnorm(mat[2,1], 0, 1)
third=rnorm(mat[3,1], 4, 2)
sample = c(rpois(n = mat[1,1],1), rpois(n = mat[2,1],50), rpois(n = mat[3,1],200))

densityplot = ggplot()+ aes(sample) + geom_density() + ggtitle("Mixture density of Y")
densityplot + ylab(expression("f"[ Y ]*"(y)"))  + xlab("Y") + theme(axis.text=element_text(size=14),axis.title=element_text(size=18), plot.title=element_text(size=20))
geom_segment(aes(x=-10, y=0, xend=-10, yend=max(first)))