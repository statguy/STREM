library(WTC)

response <- "canis.lupus"
response <- "lynx.lynx"
context <- Context$new(resultDataDirectory=wd.data.results, processedDataDirectory=wd.data.processed, rawDataDirectory=wd.data.raw, scratchDirectory=wd.scratch, figuresDirectory=wd.figures)
study <- FinlandWTCStudy$new(context=context, response=response)

populationSize <- study$getPopulationSize(withHabitatWeights=TRUE, withDistanceWeights=TRUE, saveDensityPlots=TRUE)
populationSize$loadValidationData()
populationSize$plotPopulationSize()

# Diagnostics

estimates <- study$collectEstimates()
estimates$getPredictedIntersections() # TODO: doesn't work because distances have not been eliminated

populationDensity1 <- study$getPopulationDensity(withHabitatWeights=FALSE, withDistanceWeights=FALSE, saveDensityPlots=FALSE)
populationSize1 <- populationDensity1$mean$integrate(volume=FinlandPopulationSize$new(study=study))
populationDensity2 <- study$getPopulationDensity(withHabitatWeights=TRUE, withDistanceWeights=FALSE, saveDensityPlots=FALSE)
populationSize2 <- populationDensity2$mean$integrate(volume=FinlandPopulationSize$new(study=study))
populationDensity3 <- study$getPopulationDensity(withHabitatWeights=TRUE, withDistanceWeights=TRUE, saveDensityPlots=FALSE)
populationSize3 <- populationDensity3$mean$integrate(volume=FinlandPopulationSize$new(study=study))
populationSize3$loadValidationData()

library(reshape2)
library(ggplot2)
x <- data.frame(Year=populationSize1$sizeData$Year,
                mean_distance=populationSize1$sizeData$Estimated,
                habitat_weights_mean_distance=populationSize2$sizeData$Estimated,
                habitat_weights_predicted_distances=populationSize3$sizeData$Estimated,
                adjusted=populationSize3$sizeData$Estimated,
                validation=populationSize3$sizeData$Observed)
x <- melt(x, id.vars="Year")
ggplot(x, aes(x=Year, y=value, group=variable, color=variable)) + geom_line()
