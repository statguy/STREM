# Run test:
# ./parallel_r.py -t 1:50 -n 60 -l 10.0 -b ~/tmp/blacklist.txt -v ~/git/Winter-Track-Counts/inst/simulation/estimate.R test A SmoothModel-nbinomial-matern-ar1
# Run full:
# ./parallel_r.py -t 1:50 -n 70 -l 10.0 -b ~/tmp/blacklist.txt -v ~/git/Winter-Track-Counts/inst/simulation/estimate.R notest A SmoothModel-nbinomial-matern-ar1
# ./parallel_r.py -t 1:50 -n 70 -l 10.0 -b ~/tmp/blacklist.txt -v ~/git/Winter-Track-Counts/inst/simulation/estimate.R notest A SmoothModel-nbinomial-ar1
#
# library(devtools); install_github("statguy/Winter-Track-Counts")

estimateSpatioTemporal <- function(scenario, modelName, iteration, isTest=FALSE, quick=FALSE) {
  mss <- getMSS(scenario=scenario, isTest=isTest)
  study <- mss$study
  
  if (quick) {
    stop("TODO")
    intersections <- study$loadIntersections(iteration=iteration)
    model <- SimulatedSmoothModelSpatioTemporal$new(study=study, iteration=iteration)
    meshParams <- list(maxEdge=c(.1e6, .2e6), cutOff=.05e6, coordsScale=1e-6)
    modelParams <- list(family="nbinomial", offsetScale=1000^2, meshParams=meshParams, timeModel="ar1")
    model$setup(intersections=intersections, params=modelParams)
    
    model$estimate()
    model$collectEstimates()
    model$collectHyperparameters()
    model$plotTemporal()
    
    populationSize <- model$getPopulationSize(withHabitatWeights=FALSE); populationSize
    populationSize$plotPopulationSize()
  }
  else {
    if (modelName == "SmoothModel-nbinomial-matern-ar1") {
      model <- SimulatedSmoothModelSpatioTemporal(study=study, iteration=iteration)
      modelParams <- list(family="nbinomial", offsetScale=1000^2, meshParams=study$studyArea$getMesh(), timeModel="ar1")      
    }
    else if (modelName == "SmoothModel-nbinomial-ar1") {
      model <- SimulatedSmoothModelTemporal(study=study, iteration=iteration)
      modelParams <- list(family="nbinomial", offsetScale=1000^2, timeModel="ar1")
    }
    else stop("Unknown model.")
    study$estimate(model=model, params=modelParams)
  }
}

library(parallel)
library(doMC)
registerDoMC(cores=detectCores())
library(WTC)
source("~/git/Winter-Track-Counts/setup/WTC-Boot.R")

parseArguments()
modelName <- extraArgs[1]
estimateSpatioTemporal(scenario=scenario, modelName=modelName, iteration=as.integer(task_id), isTest=isTest)
