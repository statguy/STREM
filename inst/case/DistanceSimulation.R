library(CNPCluster)
library(WTC)

cnpClusterStartLocal()

context <- Context$new(resultDataDirectory=wd.data.results, processedDataDirectory=wd.data.processed, rawDataDirectory=wd.data.raw, scratchDirectory=wd.scratch, figuresDirectory=wd.figures)
mssIntensive <- MovementSimulationScenarioIntensive$new(stepIntervalHours=4)$newInstance(context=context)
tracks <- mssIntensive$simulate(save=TRUE)

study <- mssIntensive$study
tracks <- study$loadTracksCollection()
surveyRoutes <- FinlandRandomWTCSurveyRoutes$new(study=study)$newInstance(800)
tracks$getDistances(surveyRoutes)
intersections <- tracks$findIntersections(surveyRoutes, dimension=1, save=TRUE)

intersections <- study$loadIntersectionsCollection()
#meshParams <- list(maxEdge=c(.05e6, .15e6), cutOff=.02e6, coordsScale=1e-6)
meshParams <- list(maxEdge=c(.05e6, .15e6), cutOff=.05e6, coordsScale=1e-6)
model <- FinlandSmoothModel$new(study=study)$setup(intersections=intersections$intersectionsList[[1]], meshParams=meshParams)
model$plotMesh(surveyRoutes=surveyRoutes)
model$estimate(save=TRUE)

model <- FinlandSmoothModel$new(study=study)
model$loadResult()
model$collectResults(quick=TRUE)

cnpClusterStopLocal()
