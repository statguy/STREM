library(parallel)
library(doMC)
registerDoMC(cores=detectCores())
library(WTC)
source("~/git/Winter-Track-Counts/setup/WTC-Boot.R")
library(plyr)

scenarios <- c("A","B","C","D","E","F")
allScenarios <- c(scenarios, paste0(scenarios, "combined"), paste0(scenarios, "10days"))

l_ply(allScenarios, function(scenario)) {

  mss <- getMSS(scenario=scenario)
  study <- mss$study
  
  tracks <- SimulatedTracks(study=study)
  iterations <- tracks$getTracksFileIterations()
  for (iteration in iterations) {
    tracks <- study$loadTracks(iteration=iteration, addColumns=FALSE)
    truePopulationSize <- tracks$truePopulationSize
    fileName <- study$context$getLongFileName(dir=tracks$getTracksDirectory(), name="TruePopulationSize", response=study$response, region=study$studyArea$region, tag=iteration)
    save(truePopulationSize, file=fileName)
  }

}, .parallel=T)
