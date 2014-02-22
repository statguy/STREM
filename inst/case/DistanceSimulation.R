library(devtools)
install_github("R-Cluster", "statguy")
library(devtools)
install_github("Winter-Track-Counts", "statguy")
library(devtools)
source_gist("b7507c36efada51bbda5") # TODO: remove

library(CNPCluster)
library(WTC)

cnpClusterStartLocal()

context <- Context$new(resultDataDirectory=wd.data.results, processedDataDirectory=wd.data.processed, rawDataDirectory=wd.data.raw, scratchDirectory=wd.scratch, figuresDirectory=wd.figures)
mssIntensive <- MovementSimulationScenarioIntensive$new()$newInstance(context=context)
tracks <- mssIntensive$simulate(save=TRUE)

study <- mssIntensive$study
tracks <- study$loadTracksCollection()

sampleIntervals <- MovementSampleIntervals$new(study=study)
thinnedTracks <- sampleIntervals$getThinnedTracksSampleIntervals(tracks$getTracks(1))
thinnedTracks$determineDistances()

surveyRoutes <- FinlandRandomWTCSurveyRoutes$new(study=study)$newInstance(800)
intersections <- thinnedTracks$findIntersections(surveyRoutes, dimension=1, save=FALSE)
# TODO: combine intersections here

#intersections <- study$loadIntersectionsCollection()
meshParams <- list(maxEdge=c(.05e6, .15e6), cutOff=.02e6, coordsScale=1e-6)
#meshParams <- list(maxEdge=c(.05e6, .15e6), cutOff=.05e6, coordsScale=1e-6)
models <- intersections$estimate(meshParams=meshParams, save=FALSE)

#models <- SmoothModelCollection$new(study=study, directory=study$context$scratchDirectory)$loadModels()
model <- models$getModel(1)
model$plotMesh()
model$collectResults(quick=TRUE)

habitatWeights <- HabitatWeights$new(study=study)
projectionRaster <- habitatWeights$getWeightsRaster(aggregationScale=100, save=T)
populationDensity <- model$getPopulationDensity(projectionRaster=projectionRaster, maskPolygon=NULL, getSD=FALSE)
populationSize <- populationDensity$mean$integrate(weights=1)
populationSize$sizeData


cnpClusterStopLocal()
