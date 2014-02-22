library(CNPCluster)
library(WTC)

cnpClusterStartLocal()

library(devtools)
source_gist("b7507c36efada51bbda5")
install_github("Winter-Track-Counts", username="statguy") # TODO: remove
context <- Context$new(resultDataDirectory=wd.data.results, processedDataDirectory=wd.data.processed, rawDataDirectory=wd.data.raw, scratchDirectory=wd.scratch, figuresDirectory=wd.figures)
nAgents <- as.integer(2)
mssIntensive <- MovementSimulationScenarioIntensive$new(stepIntervalHours=4, nAgents=nAgents)$newInstance(context=context)
tracks <- mssIntensive$simulate(save=TRUE)

study <- mssIntensive$study
tracks <- study$loadTracksCollection()

sampleIntervals <- MovementSampleIntervals$new(study=study)
thinnedTracks <- sampleIntervals$getThinnedTracksSampleIntervals(tracks$getTracks(1))
thinnedTracks$determineDistances()

surveyRoutes <- FinlandRandomWTCSurveyRoutes$new(study=study)$newInstance(800)
intersections <- thinnedTracks$findIntersections(surveyRoutes, dimension=1, save=FALSE)

#intersections <- study$loadIntersectionsCollection()
#meshParams <- list(maxEdge=c(.05e6, .15e6), cutOff=.02e6, coordsScale=1e-6)
meshParams <- list(maxEdge=c(.05e6, .15e6), cutOff=.05e6, coordsScale=1e-6)
models <- intersections$estimate(meshParams=meshParams, save=FALSE)

#models <- SmoothModelCollection$new(study=study, directory=study$context$scratchDirectory)$loadModels()
#models$getModel(1)$collectResults(quick=TRUE)



cnpClusterStopLocal()
