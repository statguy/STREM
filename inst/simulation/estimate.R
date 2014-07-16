# Run test:
# ./parallel_r.py -t 1:3 -n 3 -l 10.0 -b ~/tmp/blacklist.txt -v ~/git/Winter-Track-Counts/inst/simulation/estimate.R test A
# Run full:
# ./parallel_r.py -t 1:50 -n 70 -l 10.0 -b ~/tmp/blacklist.txt -v ~/git/Winter-Track-Counts/inst/simulation/estimate.R notest A

# library(devtools); install_github("statguy/Winter-Track-Counts")


estimate <- function(study, iteration, test) {  
  if (test=="test") {
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
    meshParams <- list(maxEdge=c(.1e6, .2e6), cutOff=.05e6, coordsScale=1e-6)
    modelParams <- list(family="nbinomial", offsetScale=1000^2, meshParams=meshParams, timeModel="ar1")
    study$estimate(model=model, params=modelParams)
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
mss <- {
  if (scenario == "A") MovementSimulationScenarioA$new()$setup(context=context)
  else if (scenario == "B") MovementSimulationScenarioB$new()$setup(context=context)
  else if (scenario == "C") MovementSimulationScenarioC$new()$setup(context=context)
  else if (scenario == "D") MovementSimulationScenarioD$new()$setup(context=context)
  else if (scenario == "E") MovementSimulationScenarioE$new()$setup(context=context)
  else if (scenario == "F") MovementSimulationScenarioF$new()$setup(context=context)
  else stop("unsupported")
}

study <- mss$study
estimate(study=study, iteration=as.integer(task_id), test)
