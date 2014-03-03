library(WTC)

response <- "canis.lupus"
context <- Context$new(resultDataDirectory=wd.data.results, processedDataDirectory=wd.data.processed, rawDataDirectory=wd.data.raw, scratchDirectory=wd.scratch, figuresDirectory=wd.figures)
study <- FinlandWTCStudy$new(context=context, response=response)

distances <- study$predictDistances()
populationSize <- study$getPopulationSize(weights=subset(distances, Variable=="Predicted", select="Value"))
populationSize$plotPopulationSize()
