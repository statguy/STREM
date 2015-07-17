library(parallel)
library(doMC)
registerDoMC(cores=round(detectCores()))
source("~/git/Winter-Track-Counts/setup/WTC-Boot.R")
library(WTC)
library(ggplot2)

plotMovementHabitatDistribution <- function(response) {
  context <- Context$new(resultDataDirectory=wd.data.results, processedDataDirectory=wd.data.processed, rawDataDirectory=wd.data.raw, scratchDirectory=wd.scratch, figuresDirectory=wd.figures)
  study <- FinlandWTCStudy$new(context=context)
  study$response <- response
  tracks <- study$loadTracks()
  habitatWeights <- CORINEHabitatWeights$new(study=study)
  xy <- ld(tracks$tracks)[,c("x","y")]
  habitatValues <- raster::extract(study$studyArea$habitat, xy)
  classifiedHabitatValues <- habitatWeights$weights$type[habitatValues]
  classifiedHabitatValues <- factor(classifiedHabitatValues, levels=habitatWeights$habitatTypes, labels=names(habitatWeights$habitatTypes))
  ggplot(data.frame(classifiedHabitatValues), aes(x=classifiedHabitatValues)) + geom_bar()
}

plotMovementHabitatDistribution("canis.lupus")
plotMovementHabitatDistribution("lynx.lynx")
plotMovementHabitatDistribution("rangifer.tarandus.fennicus")
