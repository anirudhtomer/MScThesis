# MScThesis
## How to read/run the code for simulation study

- Go to folder src\dic-ppc-bf-randslope
- Load the bhtge_slope.Rproj file in R.
- The file simulationHeterogeneity.R is the primary file and should be used to run all MCMC simulations.
- The JAGS related files such as model definitions and configuration of JAGS are present in src\common folder.
- The model definition is in called createModelRandSlope.R. There are separate model definitions for model with no component in mixture and models with 2 or more components. There are other files starting with the name createModel....They are used for calculation of marginal likelihood using Chib's approximation.
- fitModelRandSlope.R has the code for initialization and monitoring of JAGS. The fitModelRandSlopeConditional.R file contains a similar code but for calculation of marginal likelihood using Chib's approximation.
- The simulation study data sets can be tweaked using generateDataRandSlope.R
- The postprocessing and extraction of parameters from chain is in extractFuncRandSlope.R and in extractFuncRandSlopeConditional.R.
- The code for DIC calculation, PPC calculation and Bayes factor calculation is in src\dic-ppc-bf-randslope folder. The names of files are DIC_functions.R, ppc.R and bf.R respectively.


