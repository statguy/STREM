# library(devtools); install_github("statguy/Winter-Track-Counts")

library(parallel)
library(doMC)
registerDoMC(cores=detectCores())
library(WTC)
source("~/git/Winter-Track-Counts/setup/WTC-Boot.R")

responses <- c("canis.lupus", "lynx.lynx", "rangifer.tarandus.fennicus")
modelNames <- c("FMPModel", "SmoothModel-nbinomial-ar1", "SmoothModel-nbinomial-rw2", "SmoothModel-nbinomial-matern-ar1")

for (response in responses) {
  for (modelName in modelNames) {    
    context <- Context$new(resultDataDirectory=wd.data.results, processedDataDirectory=wd.data.processed, rawDataDirectory=wd.data.raw, scratchDirectory=wd.scratch, figuresDirectory=wd.figures)
    study <- FinlandWTCStudy$new(context=context, response=response, distanceCovariatesModel=~populationDensity+rrday+snow+tday-1, trackSampleInterval=2)
    populationSize <- study$getPopulationSize(modelName=modelName)
    
    if (F) {
      habitatWeightsRaster <- study$loadHabitatWeightsRaster()
      model <- study$getModel(modelName=modelName)  
      model$offsetScale <- 1000^2 # quickfix, remove when not needed anymore
      model$loadEstimates()
      model$collectEstimates()
      populationDensity <- model$getPopulationDensity(templateRaster=habitatWeightsRaster, getSD=FALSE)
      populationSize <- FinlandPopulationSize(study=study, modelName=modelName)$getPopulationSize(populationDensity=populationDensity$mean)
    }
    
    populationSize$plotPopulationSize()
    print(populationSize)
  }
}