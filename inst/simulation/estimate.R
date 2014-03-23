# Run test:
# ./parallel_r.py -t 1:3 -n 3 -l 10.0 -b ~/tmp/blacklist.txt -v ~/git/Winter-Track-Counts/inst/simulation/estimate.R test A
# Run full:
# ./parallel_r.py -t 1:50 -n 60 -l 10.0 -b ~/tmp/blacklist.txt -v ~/git/Winter-Track-Counts/inst/simulation/estimate.R notest A

# library(devtools); install_github("statguy/Winter-Track-Counts")


estimate <- function(study, iteration, test) {  
  if (test) {
    intersections <- study$loadIntersections(iteration=iteration)
    tracks <- study$loadTracks(iteration=iteration)
    intersections$intersections$distance <- tracks$getMeanDistance()
    model <- SimulatedSmoothModel$new(study=study, iteration=iteration)
    meshParams <- list(maxEdge=c(.2e6, .3e6), cutOff=.2e6, coordsScale=1e-6)
    model$setup(intersections=intersections, meshParams=meshParams)
    model$estimate()
    model$collectEstimates()
    
    #sum(model$node$mean[,1] * (inla.mesh.fem(model$mesh, order=1)$c1 %*% model$node$mean[,1]))
  }
  else {
    intersections <- study$loadIntersections(iteration=iteration)
    #meshParams <- list(maxEdge=c(.05e6, .15e6), cutOff=.02e6, coordsScale=1e-6)
    meshParams <- list(maxEdge=c(.1e6, .2e6), cutOff=.05e6, coordsScale=1e-6)
    model <- intersections$estimate(meshParams=meshParams)
    model$saveEstimates()
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
