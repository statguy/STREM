# Run test:
# ./parallel_r.py -t 1:3 -n 3 -l 10.0 -b ~/tmp/blacklist.txt -v ~/git/Winter-Track-Counts/inst/simulation/estimate.R test A
# Spatial test:
# ./parallel_r.py -t 1:50 -n 70 -l 10.0 -b ~/tmp/blacklist.txt -v ~/git/Winter-Track-Counts/inst/simulation/estimate.R spatial A
# Run full:
# ./parallel_r.py -t 1:50 -n 70 -l 10.0 -b ~/tmp/blacklist.txt -v ~/git/Winter-Track-Counts/inst/simulation/estimate.R notest A

# library(devtools); install_github("statguy/Winter-Track-Counts")


estimate <- function(study, iteration, test) {  
  if (test=="test") {
    intersections <- study$loadIntersections(iteration=iteration)
    model <- SimulatedSmoothModel$new(study=study, iteration=iteration)
    #model <- SimulatedSmoothModelNoOffset$new(study=study, iteration=iteration)
    #offsets <- SimulatedSmoothModel$new(study=study, iteration=iteration)
    meshParams <- list(maxEdge=c(.2e6, .3e6), cutOff=.2e6, coordsScale=1e-6)
    #meshParams <- list(maxEdge=c(.07e6, .2e6), cutOff=.01e6, coordsScale=1e-6)
    interceptPriorParams <- list(mean=200, sd=199)
    rhoPriorParams <- list(initial=0, param=c(0, 0.5))
    model$setup(intersections=intersections, meshParams=meshParams, useCovariates=FALSE)
    model$setupInterceptPrior(interceptPriorParams)
    model$setupRhoPrior(rhoPriorParams)
    c(mean(model$getObservedOffset()),sd(model$getObservedOffset()))
    c(mean(model$getPredictedOffset()),sd(model$getPredictedOffset()))
    model$estimate()
    model$collectEstimates()
    #model$collectEstimates(weightsAtSurveyRoutes=offsets$getObservedOffset(), weightsAtNodes=offsets$getNodeOffset())
    model$collectHyperparameters()
    model$plotTemporal()
    
    populationSize <- model$getPopulationSize(withHabitatWeights=FALSE); populationSize
    populationSize$plotPopulationSize()
    
    #sum(model$node$mean[,1] * (inla.mesh.fem(model$mesh, order=1)$c1 %*% model$node$mean[,1]))
  }
  else if (test=="spatial") {
    message("Spatial test...")
    meshParams <- list(maxEdge=c(.1e6, .2e6), cutOff=.05e6, coordsScale=1e-6)
    interceptPriorParams <- list(mean=200, sd=199)
    study <- SimulationStudySubset$new(response="A", years=as.integer(2001:2001))$newInstance(context=context)
    intersections <- study$loadIntersections(iteration=iteration)
    model <- SimulatedSmoothModel$new(study=study, iteration=iteration)
    model$setup(intersections=intersections, meshParams=meshParams, useCovariates=FALSE)
    model$setupInterceptPrior(interceptPriorParams)
    model$estimate()
    model$saveEstimates()
  }
  else {
    meshParams <- list(maxEdge=c(.1e6, .2e6), cutOff=.05e6, coordsScale=1e-6)
    interceptPriorParams <- list(mean=200, sd=199)
    study$estimate(iteration=iteration, meshParams=meshParams, interceptPriorParams=interceptPriorParams)
  }
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 3) stop("Invalid arguments.")
test <- args[1]
scenario <- args[2]
task_id <- args[length(args)]
message("Arguments provided:")
print(args)

library(parallel)
library(doMC)
registerDoMC(cores=detectCores())
library(WTC)
source("~/git/Winter-Track-Counts/setup/WTC-Boot.R")

context <- Context$new(resultDataDirectory=wd.data.results, processedDataDirectory=wd.data.processed, rawDataDirectory=wd.data.raw, scratchDirectory=wd.scratch, figuresDirectory=wd.figures)
mss <- if (scenario == "A") MovementSimulationScenarioA$new()$newInstance(context=context) else stop("Unknown scenario ", scenario)
study <- mss$study
estimate(study=study, iteration=as.integer(task_id), test)
