# library(devtools); install_github("statguy/Winter-Track-Counts")

library(parallel)
library(doMC)
registerDoMC(cores=detectCores())
library(WTC)
source("~/git/Winter-Track-Counts/setup/WTC-Boot.R")

task_id <- 5

nSurveyRoutes <- as.integer(100)
nAgents <- as.integer(20)
nIterations <- as.integer(100)
nYears <- as.integer(1)
iteration <- as.integer(1)
nDays <- as.integer(59)
meshParams <- list(coordsScale=1e-6, maxEdge=c(.01e6, .02e6), cutOff=.007e6)
modelParams <- list(family="nbinomial", offsetScale=1000^2, meshParams=meshParams, timeModel="ar1")
BCRWCorrelationBiasTradeoff <- 0.7

context <- Context$new(resultDataDirectory=wd.data.results, processedDataDirectory=wd.data.processed, rawDataDirectory=wd.data.raw, scratchDirectory=wd.scratch, figuresDirectory=wd.figures)
{
mss <- if (task_id == 1) MovementSimulationScenarioA$new(nAgents=nAgents, years=nYears, days=nDays)$newInstance(context=context, isTest=T)
else if (task_id == 2) MovementSimulationScenarioB$new(nAgents=nAgents, years=nYears, days=nDays, BCRWCorrelationBiasTradeoff=BCRWCorrelationBiasTradeoff)$newInstance(context=context, isTest=T)

else if (task_id == 4) MovementSimulationScenarioD$new(nAgents=nAgents, years=nYears, days=nDays)$newInstance(context=context, isTest=T)
else if (task_id == 5) MovementSimulationScenarioE$new(nAgents=nAgents, years=nYears, days=nDays, nSurveyRoutes=nSurveyRoutes)$newInstance(context=context, isTest=T)

}

study <- mss$study
message("Study area = ", study$studyArea$boundary@polygons[[1]]@area / 1000^2, " km^2")
habitatWeights <- CORINEHabitatWeights$new(study=study)
size <- c()

for (iteration in 1:nIterations) {
  #mss$debug <- TRUE
  tracks <- mss$simulate(iteration=iteration, save=F)
  surveyRoutes <- mss$getSurveyRoutes(nSurveyRoutes=nSurveyRoutes)
  #tracks$plotTracks(surveyRoutes=surveyRoutes, habitat=T)
  #plot(mss$initialPopulation$locations, add=T)
  
  intersections <- SimulatedIntersections$new(study=study, iteration=iteration)
  intersections$findIntersections(tracks, surveyRoutes,  dimension=1)
  
  model <- SimulatedSmoothModelSpatioTemporal$new(study=study, iteration=iteration)
  model$setup(intersections=intersections, params=modelParams)
  #plot(model$getUnscaledMesh()); plot(study$studyArea$boundary, add=T, border="blue")
  model$estimate()
  model$collectEstimates()
  
  habitatSelection <- tracks$getHabitatPreferences(habitatWeightsTemplate=habitatWeights, nSamples=30, save=FALSE)
  habitatWeights$setHabitatSelectionWeights(habitatSelection)
  habitatWeightsRaster <- habitatWeights$getWeightsRaster(save=FALSE)
  #habitatWeights; plot(habitatWeightsRaster)
  populationDensity <- model$getPopulationDensity(getSD=FALSE)
  
  populationSize <- populationDensity$mean$integrate(volume=SimulationPopulationSize$new(study=study, iteration=iteration))
  populationDensity$mean$weight(habitatWeightsRaster)
  populationSizeWeighted <- populationDensity$mean$integrate(volume=SimulationPopulationSize$new(study=study, iteration=iteration))
  
  habitatWeights; populationSize; populationSizeWeighted
  size <- rbind(size, data.frame(obs=sum(model$data$intersections), fitted=sum(model$data$fittedMean * model$getObservedOffset()), sizew=populationSizeWeighted$sizeData$Estimated, size=populationSize$sizeData$Estimated, iteration=iteration))
  
  
  #populationSize <- model$getPopulationSize(tracks=tracks, withHabitatWeights=mss$hasHabitatWeights())
  #x <- populationSize$sizeData$Estimated
  #model$switchToMesh()
  #populationSize <- model$getPopulationSize(tracks=tracks, withHabitatWeights=mss$hasHabitatWeights())
  
  #size <- rbind(size, data.frame(obs=sum(model$data$intersections), fitted=sum(model$data$fittedMean * model$getObservedOffset()), sizeobs=x, sizenode=populationSize$sizeData$Estimated, iteration=iteration))
}

colMeans(size); colSDs(size)
