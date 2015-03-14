# Run test:
# ./parallel_r.py -t 1:50 -n 60 -l 10.0 -b ~/tmp/blacklist.txt -v ~/git/Winter-Track-Counts/inst/simulation/estimate.R test A SmoothModel-nbinomial-matern-ar1
# Run full:
# ./parallel_r.py -t 1:50 -n 70 -l 10.0 -b ~/tmp/blacklist.txt -v ~/git/Winter-Track-Counts/inst/simulation/estimate.R notest A SmoothModel-nbinomial-matern-ar1
# ./parallel_r.py -t 1:50 -n 70 -l 10.0 -b ~/tmp/blacklist.txt -v ~/git/Winter-Track-Counts/inst/simulation/estimate.R notest A SmoothModel-nbinomial-ar1
# ./parallel_r.py -t 1:50 -n 70 -l 10.0 -b ~/tmp/blacklist.txt -v ~/git/Winter-Track-Counts/inst/simulation/estimate.R notest A FMPModel
# ./parallel_r.py -t 1:5 -n 7 -l 10.0 -b ~/tmp/blacklist.txt -v ~/git/Winter-Track-Counts/inst/simulation/estimate.R notest Acombined SmoothModel-nbinomial-matern-ar1
# ./parallel_r.py -t 1:5 -n 7 -l 10.0 -b ~/tmp/blacklist.txt -v ~/git/Winter-Track-Counts/inst/simulation/estimate.R notest Acombined SmoothModel-nbinomial-ar1
# ./parallel_r.py -t 1:5 -n 7 -l 10.0 -b ~/tmp/blacklist.txt -v ~/git/Winter-Track-Counts/inst/simulation/estimate.R notest Acombined FMPModel
#
# library(devtools); install_github("statguy/Winter-Track-Counts")

estimateSpatioTemporal <- function(scenario, modelName, iteration, isTest=FALSE, quick=FALSE) {
  mss <- getMSS(scenario=scenario, isTest=isTest)
  study <- mss$study
  
  if (quick) {
    stop("TODO")
    intersections <- study$loadIntersections(iteration=iteration)
    model <- SimulatedSmoothModelSpatioTemporal(study=study, iteration=iteration)
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
    #model <- study$getModel(modelName=modelName)
    #modelParams <- study$getModelParams(modelName=modelName)
    tag <- NULL
    
    if (modelName == "SmoothModel-nbinomial-matern-ar1") {
      model <- SimulatedSmoothModelSpatioTemporal(study=study, iteration=iteration)
      modelParams <- list(family="nbinomial", offsetScale=1000^2, meshParams=study$studyArea$getMesh(), timeModel="ar1")
    }
    else if (modelName == "SmoothModel-nbinomial-ar1") {
      model <- SimulatedSmoothModelTemporal(study=study, iteration=iteration)
      formula <- response ~ 1 + f(year, model="ar1")
      modelParams <- list(family="nbinomial", offsetScale=1000^2, model=formula, timeModel="ar1")
    }
    else if (modelName == "SmoothModelMean-nbinomial-ar1") {
      model <- SimulatedSmoothModelMeanTemporal(study=study, iteration=iteration)
      formula <- response ~ 1 + f(year, model="ar1")
      modelParams <- list(family="nbinomial", offsetScale=1000^2, model=formula, timeModel="ar1")
    }
    else if (modelName == "SmoothModelMean-nbinomial-ar1-priors1") {
      model <- SimulatedSmoothModelMeanTemporal(study=study, iteration=iteration)
      #precPrior <- model$setupPrecisionPrior(priorParams=list(shape=1, scale=5e-5, initial=4))
      #rhoPrior <- model$setupTemporalPrior(priorParams=list(mean=0, sd=0.15, initial=2))
      #precPrior <- model$setupPrecisionPrior(priorParams=list(mean=1, sd=1, initial=4))
      #rhoPrior <- model$setupTemporalPrior(priorParams=list(shape=0, rate=1, initial=2))
      #formula <- response ~ 1 + f(year, model="ar1", hyper=list(theta1=precPrior, theta2=rhoPrior))
      formula <- response ~ 1 + f(year, model="ar1", hyper=list(theta1=list(1, 1), theta2=list(0, 1)))
      modelParams <- list(family="nbinomial", offsetScale=1000^2, model=formula, timeModel="ar1")
      tag <- "priors1"
    }
    else if (modelName == "FMPModel") {
      model <- SimulatedFMPModel(study=study, iteration=iteration)
      modelParams <- NULL
    }
    else stop("Unknown model.")
    study$estimate(model=model, params=modelParams, tag=tag)
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
#estimateSpatioTemporal(scenario="A", modelName="SmoothModelMean-nbinomial-ar1", iteration=as.integer(1), isTest=FALSE)
