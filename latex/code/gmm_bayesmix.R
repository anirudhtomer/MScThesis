install.packages("bayesmix")
library(bayesmix)

model <- BMMmodel(sample, k = 3, initialValues = list(S0 = 2),
                  priors = list(kind = "independence",
                  parameter = "priorsFish", hierarchical = "tau"))  
control <- JAGScontrol(variables = c("mu", "tau", "eta", "S"),
                      burn.in = 1000, n.iter = 10000, seed = 10)
z <- JAGSrun(sample, model = model, control = control)
BMMdiag(z)
