library(devtools)
install_github("statguy/WTC")
install_github("ropengov/gisfin")
install_github("ropengov/fmi")

library(parallel)
library(doMC)
registerDoMC(cores=round(detectCores()))
source("~/git/Winter-Track-Counts/setup/WTC-Boot.R")
library(WTC)

context <- Context$new(resultDataDirectory=wd.data.results, processedDataDirectory=wd.data.processed, rawDataDirectory=wd.data.raw, scratchDirectory=wd.scratch, figuresDirectory=wd.figures)
study <- FinlandWTCStudy$new(context=context)

if (F) {
  study$response <- "canis.lupus"
  
  ###
  
  tracks <- study$loadTracks()
  habitatWeights <- CORINEHabitatWeights$new(study=study)
  habitatSelection <- tracks$getHabitatPreferences(habitatWeightsTemplate=habitatWeights, nSamples=30, save=F)
  habitatWeights <- CORINEHabitatWeights$new(study=study)$setHabitatSelectionWeights(habitatSelection)
  
  ###
  
  tracks <- study$loadTracks()
  sampleIntervals <- ThinnedMovementSampleIntervals$new()
  sampleIntervals$findSampleIntervals(tracks=tracks$tracks)
  
  ###
  
  xyt <- sampleIntervals$getSampleLocations()
  
  covariates <- HumanPopulationDensityCovariates$new(study=study)
  covariates$preprocess()
  human <- covariates$extract(xyt)
  
  covariates <- WeatherCovariates$new(study=study, apiKey=fmiApiKey)
  covariates$preprocess()
  weather <- covariates$extract(xyt)
  
  sampleIntervals$associateCovariates(human, weather)
  sampleIntervals$save()
}
else {
  study$preprocess(fmiApiKey=fmiApiKey) 
}
