library(parallel)
library(doMC)
registerDoMC(cores=detectCores())
library(WTC)
source("~/git/Winter-Track-Counts/setup/WTC-Boot.R")

task_id <- 1

nSurveyRoutes <- as.integer(50)
nAgents <- as.integer(50)
nIterations <- as.integer(2)
nYears <- as.integer(1)
iteration <- as.integer(1)
nDays <- as.integer(59)
meshParams <- list(coordsScale=1e-6, maxEdge=c(.01e6, .02e6), cutOff=.005e6)
modelParams <- list(family="nbinomial", offsetScale=1000^2, meshParams=meshParams, timeModel="ar1")
BCRWCorrelationBiasTradeoff <- 0.7

context <- Context$new(resultDataDirectory=wd.data.results, processedDataDirectory=wd.data.processed, rawDataDirectory=wd.data.raw, scratchDirectory=wd.scratch, figuresDirectory=wd.figures)
mss <- if (task_id == 1) MovementSimulationScenarioB$new(nAgents=nAgents, years=nYears, days=nDays)$newInstance(context=context, isTest=T)
else if (task_id == 2) MovementSimulationScenarioB$new(nAgents=nAgents, years=nYears, days=nDays, BCRWCorrelationBiasTradeoff=BCRWCorrelationBiasTradeoff)$newInstance(context=context, isTest=T)
study <- mss$study
surveyRoutes <- FinlandRandomWTCSurveyRoutes$new(study=study)$randomizeSurveyRoutes(nSurveyRoutes=nSurveyRoutes)

populationSizes <- c()

for (iteration in 1:nIterations) {
  tracks <- mss$simulate(iteration=iteration, save=F)
  #tracks$plotTracks(surveyRoutes=surveyRoutes)

  intersections <- SimulatedIntersections$new(study=study, iteration=iteration)
  intersections$findIntersections(tracks, surveyRoutes,  dimension=1)
  
  model <- SimulatedSmoothModelSpatioTemporal$new(study=study, iteration=iteration)
  model$setup(intersections=intersections, params=modelParams)
  #plot(model$mesh)
  model$estimate()
  model$collectEstimates()
  
  model$switchToMesh()
  populationSize <- model$getPopulationSize(tracks=tracks, withHabitatWeights=FALSE)
  x <- populationSize$sizeData$Estimated
  message("Iteration ", iteration, " / ", nIterations, ", n = ", x)
  populationSizes <- c(populationSizes, x)
}

mean(populationSizes); sd(populationSizes)
