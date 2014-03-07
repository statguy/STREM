library(WTC)

response <- "lynx.lynx"
context <- Context$new(resultDataDirectory=wd.data.results, processedDataDirectory=wd.data.processed, rawDataDirectory=wd.data.raw, scratchDirectory=wd.scratch, figuresDirectory=wd.figures)
study <- FinlandWTCStudy$new(context=context, response=response)

distances <- study$predictDistances()
populationSize <- study$getPopulationSize(volome=FinlandPopulationSize$new(study=study), weights=subset(distances, Variable=="Predicted", select="Value", drop=TRUE))
populatuonSize$loadValidationData()
populationSize$plotPopulationSize()
