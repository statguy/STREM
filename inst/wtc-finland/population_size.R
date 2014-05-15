library(parallel)
library(doMC)
registerDoMC(cores=detectCores())
library(WTC)
source("~/git/Winter-Track-Counts/setup/WTC-Boot.R")


response <- "canis.lupus"
response <- "lynx.lynx"
response <- "rangifer.tarandus.fennicus"
context <- Context$new(resultDataDirectory=wd.data.results, processedDataDirectory=wd.data.processed, rawDataDirectory=wd.data.raw, scratchDirectory=wd.scratch, figuresDirectory=wd.figures)
study <- FinlandWTCStudy$new(context=context, response=response, distanceCovariatesModel=~populationDensity+rrday+snow+tday-1, trackSampleInterval=2)

#populationSize <- study$getPopulationSize(withHabitatWeights=TRUE, saveDensityPlots=TRUE)
populationSize <- study$getPopulationSize(withHabitatWeights=TRUE, saveDensityPlots=FALSE)
populationSize$loadValidationData()
populationSize$plotPopulationSize()
populationSize



if (F) {
  study <- FinlandWTCStudy$new(context=context, response=response, distanceCovariatesModel=~populationDensity+rrday+snow+tday-1, trackSampleInterval=2)
  estimates <- study$loadEstimates()
  estimates$offsetScale <- 1000^2
  estimates$collectEstimates()
  summary(estimates$result)
  estimates$collectHyperparameters()
  habitatWeights <- study$loadHabitatWeightsRaster()
  populationDensity <- estimates$getPopulationDensity(templateRaster=habitatWeights, getSD=FALSE)
  populationDensity$mean$weight(habitatWeights)
  populationSize <- populationDensity$mean$integrate(volume=FinlandPopulationSize$new(study=study))
  #estimates$collectHyperparameters()
  populationSize$loadValidationData()
  populationSize
}
