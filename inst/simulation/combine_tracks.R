# library(devtools); install_github("statguy/Winter-Track-Counts")

library(parallel)
library(doMC)
registerDoMC(cores=detectCores())
library(WTC)
source("~/git/Winter-Track-Counts/setup/WTC-Boot.R")

combineTracks <- function(scenario, combineFactor) {
  mss <- getMSS(scenario=scenario)
  study <- mss$study
  tracks <- SimulatedTracks(study=study)
  iterations <- sort(tracks$getTracksFileIterations())
  end <- 1:floor(length(iterations) / combineFactor) * combineFactor
  start <- end - (combineFactor - 1)
  
  for (combineIndex in 1:length(start)) {
    iterationsIndex <- start[combineIndex]:end[combineIndex]
    message("New iteration = ", combineIndex, "/", length(start))
    message("Combine index = ", paste(iterationsIndex, collapse=" "))
    message("Combined iterations = ", paste(iterations[iterationsIndex], collapse=" "))
    
    mss <- getMSS(scenario=scenario)
    studyCombined <- mss$study
    studyCombined$response <- paste0(scenario, "combined")
    tracksCombined <- SimulatedTracks(study=studyCombined, iteration=combineIndex)
    tracksCombined$tracks <- data.frame()
    
    lastId <- 0
    for (iteration in iterations[iterationsIndex]) {
      message("Processing iteration ", iteration, "...")
      tracks <- study$loadTracks(iteration=iteration)
      tracks$tracks$id <- tracks$tracks$id + lastId
      tracks$tracks$burst <- with(tracks$tracks, paste(id, year))
      lastId <- max(tracks$tracks$id)
      tracksCombined$tracks <- rbind(tracksCombined$tracks, tracks$tracks)
    }
    
    tracksCombined$saveTracks()
  }
}

combineFactor <- 10
combineTracks(scenario="A", combineFactor=combineFactor)
combineTracks(scenario="B", combineFactor=combineFactor)
combineTracks(scenario="C", combineFactor=combineFactor)
combineTracks(scenario="D", combineFactor=combineFactor)
combineTracks(scenario="E", combineFactor=combineFactor)
combineTracks(scenario="F", combineFactor=combineFactor)
