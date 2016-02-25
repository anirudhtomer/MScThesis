install.packages("bayesmix")
library(bayesmix)

model <- BMMmodel(sample, k = 3, initialValues = list(S0 = 2),
                  priors = list(kind = "independence",
                  parameter = list("b0"=0, "B0"=10000, "nu0"=0.005, "S0"=1), hierarchical = NULL))  
control <- JAGScontrol(variables = c("mu", "tau", "eta", "S"),
                      burn.in = 1000, n.iter = 10000, seed = 10)
z <- JAGSrun(sample, model = model, control = control)
BMMdiag(z)
