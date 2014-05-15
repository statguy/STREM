library(parallel)
library(doMC)
registerDoMC(cores=detectCores())
library(WTC)
source("~/git/Winter-Track-Counts/setup/WTC-Boot.R")


response <- "canis.lupus"
response <- "lynx.lynx"
response <- "rangifer.tarandus.fennicus"
context <- Context$new(resultDataDirectory=wd.data.results, processedDataDirectory=wd.data.processed, rawDataDirectory=wd.data.raw, scratchDirectory=wd.scratch, figuresDirectory=wd.figures, verbose=TRUE)
study <- FinlandRussiaWTCStudy$new(context=context, response=response)

#populationSize <- study$getPopulationSize(withHabitatWeights=TRUE, saveDensityPlots=TRUE)
populationSize <- study$getPopulationSize()
populationSize$loadValidationData()
populationSize$plotPopulationSize()
populationSize
