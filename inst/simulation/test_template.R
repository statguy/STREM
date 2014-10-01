library(devtools); install_github("statguy/Winter-Track-Counts")

library(parallel)
library(doMC)
registerDoMC(cores=detectCores())
library(WTC)
source("~/git/Winter-Track-Counts/setup/WTC-Boot.R")

scenarios <- c("A","B","C","D","E","F")
modelNames <- c("FMPModel","SmoothModel-nbinomial-ar1","SmoothModel-nbinomial-matern-ar1")
isTest <- F
iteration <- as.integer(1)

scenario <- scenarios[6]
modelName <- modelNames[2]
mss <- getMSS(scenario=scenario, isTest=isTest)
study <- mss$study

estimates <- study$getModel(modelName=modelName, iteration=iteration)
estimates <- study$loadEstimates(estimates=estimates)
estimates$collectEstimates()
data <- estimates$data
offsetScale <- estimates$offsetScale
