library(parallel)
library(doMC)
registerDoMC(cores=detectCores())

library(WTC)
source("~/git/Winter-Track-Counts/setup/WTC-Boot.R")

n <- 2
response <- "A"
  
context <- Context$new(resultDataDirectory=wd.data.results, processedDataDirectory=wd.data.processed, rawDataDirectory=wd.data.raw, scratchDirectory=wd.scratch, figuresDirectory=wd.figures)
mss <- MovementSimulationScenarioCombined$new(n=as.integer(n))$newInstance(context=context, response=response)
mss$combine()


# Test
if (F) {
  iteration <- as.integer(2)
  study <- mss$study
  tracks <- study$loadTracks(iteration=iteration)
  intersections <- study$loadIntersections(iteration=iteration)
  meshParams <- list(maxEdge=c(.2e6, .3e6), cutOff=.2e6, coordsScale=1e-6)
  estimates <- study$estimate(meshParams=meshParams, iteration=iteration, save=FALSE)
  #populationSize <- study$getPopulationSize()

  estimates$collectEstimates()
  estimates$collectHyperparameters()
  populationSize <- estimates$getPopulationSize()
  populationSize$loadValidationData()
  populationSize
}
