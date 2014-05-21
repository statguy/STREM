library(parallel)
library(doMC)
registerDoMC(cores=detectCores())
library(WTC)
source("~/git/Winter-Track-Counts/setup/WTC-Boot.R")

responses <- c("canis.lupus", "lynx.lynx", "rangifer.tarandus.fennicus")
response <- responses[task_id]
timeModels <- c("ar1", "ar1", "rw2")

for (i in 1:length(responses)) {
  response <- responses[i]
  timeModel <- timeModels[i]
  
  context <- Context$new(resultDataDirectory=wd.data.results, processedDataDirectory=wd.data.processed, rawDataDirectory=wd.data.raw, scratchDirectory=wd.scratch, figuresDirectory=wd.figures)
  study <- FinlandWTCStudy$new(context=context, response=response, distanceCovariatesModel=~populationDensity+rrday+snow+tday-1, trackSampleInterval=2)
  model <- FinlandSmoothModelTemporal(study=study)$setModelName("nbinomial", timeModel)
  
  #populationSize <- study$getPopulationSize(model=model, withHabitatWeights=TRUE, saveDensityPlots=TRUE)
  populationSize <- study$getPopulationSize(model=model, withHabitatWeights=TRUE, saveDensityPlots=FALSE)
  populationSize$loadValidationData()
  populationSize$plotPopulationSize()
  print(populationSize)
}
