# Run test:
# ./parallel_r.py -t 1:50 -n 60 -l 10.0 -b ~/tmp/blacklist.txt -v ~/git/Winter-Track-Counts/inst/simulation/estimate.R test A
# Run full:
# ./parallel_r.py -t 1:50 -n 70 -l 10.0 -b ~/tmp/blacklist.txt -v ~/git/Winter-Track-Counts/inst/simulation/estimate.R notest A

# library(devtools); install_github("statguy/Winter-Track-Counts")


estimateSpatioTemporal <- function(scenario, iteration, isTest=FALSE, quick=FALSE) {
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
    model <- SimulatedSmoothModelSpatioTemporal(study=study, iteration=iteration)
    modelParams <- list(family="nbinomial", offsetScale=1000^2, meshParams=study$studyArea$getMesh(), timeModel="ar1")
    study$estimate(model=model, params=modelParams)
  }
}

library(parallel)
library(doMC)
registerDoMC(cores=detectCores())
library(WTC)
source("~/git/Winter-Track-Counts/setup/WTC-Boot.R")

parseArguments()

estimateSpatioTemporal(scenario=scenario, iteration=as.integer(task_id), isTest=isTest)
