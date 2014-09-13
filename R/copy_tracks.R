# library(devtools); install_github("statguy/Winter-Track-Counts")

library(parallel)
library(doMC)
registerDoMC(cores=detectCores())
library(WTC)
source("~/git/Winter-Track-Counts/setup/WTC-Boot.R")

copyTracks <- function(scenario, suffix, maxIterations) {
  mss <- getMSS(scenario=scenario)
  study <- mss$study
  tracks <- SimulatedTracks(study=study)
  iterations <- sort(tracks$getTracksFileIterations())
  iterations <- iterations[1:min(length(iterations), maxIterations)]
  message("Copying iterations ", paste(iterations, collapse=" "))
  
  for (iteration in iterations) {
    message("Iteration ", iteration, " / ", max(iterations))
    mss <- getMSS(scenario=scenario)
    studyCopy <- mss$study
    studyCopy$response <- paste0(scenario, suffix)
    tracks <- study$loadTracks(iteration=iteration)
    tracks$study <- studyCopy
    tracks$saveTracks()
  }
}

suffix <- "10days"
maxIterations<-50
copyTracks(scenario="A", suffix=suffix, maxIterations=maxIterations)
copyTracks(scenario="B", suffix=suffix, maxIterations=maxIterations)
copyTracks(scenario="C", suffix=suffix, maxIterations=maxIterations)
copyTracks(scenario="D", suffix=suffix, maxIterations=maxIterations)
copyTracks(scenario="E", suffix=suffix, maxIterations=maxIterations)
copyTracks(scenario="F", suffix=suffix, maxIterations=maxIterations)
