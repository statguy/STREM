# Run test:
# ./parallel_r.py -t 1:3 -n 3 -l 10.0 -b ~/tmp/blacklist.txt -v ~/git/Winter-Track-Counts/inst/simulation/estimate.R test A
# Run full:
# ./parallel_r.py -t 1:50 -n 60 -l 10.0 -b ~/tmp/blacklist.txt -v ~/git/Winter-Track-Counts/inst/simulation/estimate.R notest A

# library(devtools); install_github("statguy/Winter-Track-Counts")


estimate <- function(study, iteration, test) {  
  if (test) {
    intersections <- study$loadIntersections(iteration=iteration)
    model <- SimulatedSmoothModel$new(study=study, iteration=iteration)
    #model <- SimulatedSmoothModelNoOffset$new(study=study, iteration=iteration)
    meshParams <- list(maxEdge=c(.2e6, .3e6), cutOff=.2e6, coordsScale=1e-6)
    #meshParams <- list(maxEdge=c(.07e6, .2e6), cutOff=.01e6, coordsScale=1e-6)
    model$setup(intersections=intersections, meshParams=meshParams, useCovariates=FALSE)
    model$getObservedOffset(); model$getNodeOffset()
    model$estimate()
    model$collectEstimates()
    
    populationSize <- model$getPopulationSize(withHabitatWeights=FALSE); populationSize
    populationSize$plotPopulationSize()
    
    #sum(model$node$mean[,1] * (inla.mesh.fem(model$mesh, order=1)$c1 %*% model$node$mean[,1]))
  }
  else {
    meshParams <- list(maxEdge=c(.1e6, .2e6), cutOff=.05e6, coordsScale=1e-6)
    interceptPriorParams <- list(mean=200, sd=199)
    study$estimate(iteration=iteration, meshParams=meshParams, interceptPriorParams=interceptPriorParams, useCovariates=FALSE)
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
estimate(study=study, iteration=as.integer(task_id), test=test=="test")
