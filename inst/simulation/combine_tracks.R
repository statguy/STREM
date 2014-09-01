# library(devtools); install_github("statguy/Winter-Track-Counts")

library(parallel)
library(doMC)
registerDoMC(cores=detectCores())
library(WTC)
source("~/git/Winter-Track-Counts/setup/WTC-Boot.R")

combineFactor <- 10
scenario <- "A"
mss <- getMSS(scenario=scenario)
study <- mss$study

tracks <- SimulatedTracks(study=study)
iterations <- sort(tracks$getTracksFileIterations())
end <- 1:floor(length(iterations) / combineFactor) * combineFactor
start <- end - (combineFactor - 1)
for (combineIndex in 1:length(start)) {
  iterationsIndex <- start[combineIndex]:end[combineIndex]
  message("Combine index = ", paste(iterationsIndex, collapse=" "))
  message("Iterations = ", paste(iterations[iterationsIndex], collapse=" "))
  
  mss <- getMSS(scenario=scenario)
  studyCombined <- mss$study
  studyCombined$response <- paste0(scenario, "combined")
  tracksCombined <- SimulatedTracks(study=studyCombined)
  
  lastId <- 0
  for (iteration in iterations[iterationsIndex]) {
    message("Processing iteration ", iteration, "...")
    tracks <- study$loadTracks(iteration=iteration)
    tracks$tracks$id <- tracks$tracks$id + lastId
    #tracks$tracks$burst <- # TODO
    lastId <- max(tracks$tracks$id)
    tracksCombined$tracks <- rbind(tracksCombined$tracks, tracks$tracks)
    # TODO
  }
  
  tracksCombined$saveTracks()
}
