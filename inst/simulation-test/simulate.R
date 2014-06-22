# library(devtools); install_github("statguy/Winter-Track-Counts")

simulate = function(scenario, nSurveyRoutes=as.integer(50), nAgents=as.integer(20), nYears=as.integer(1),
                    nDays=as.integer(59), nIterations=as.integer(100),
                    meshParams=list(coordsScale=1e-6, maxEdge=c(.01e6, .02e6), cutOff=.007e6),
                    modelParams=list(family="nbinomial", offsetScale=1000^2, meshParams=meshParams, timeModel="ar1"),
                    CRWCorrelation=0.7,
                    BCRWCorrelationBiasTradeoff=0.7,
                    plot=FALSE) {
  library(WTC)
  source("~/git/Winter-Track-Counts/setup/WTC-Boot.R")
  
  context <- Context$new(resultDataDirectory=wd.data.results, processedDataDirectory=wd.data.processed, rawDataDirectory=wd.data.raw, scratchDirectory=wd.scratch, figuresDirectory=wd.figures)
  mss <- {
    if (scenario == "A") MovementSimulationScenarioA$new(nAgents=nAgents, years=nYears, days=nDays, CRWCorrelation=CRWCorrelation)$newInstance(context=context, isTest=T)
    else if (scenario == "B") MovementSimulationScenarioB$new(nAgents=nAgents, years=nYears, days=nDays, CRWCorrelation=CRWCorrelation, BCRWCorrelationBiasTradeoff=BCRWCorrelationBiasTradeoff)$newInstance(context=context, isTest=T)
    else if (scenario == "C") MovementSimulationScenarioC$new(nAgents=nAgents, years=nYears, days=nDays, CRWCorrelation=CRWCorrelation)$newInstance(context=context, isTest=T)
    else if (scenario == "D") MovementSimulationScenarioD$new(nAgents=nAgents, years=nYears, days=nDays, CRWCorrelation=CRWCorrelation)$newInstance(context=context, isTest=T)
    else if (scenario == "E") MovementSimulationScenarioE$new(nAgents=nAgents, years=nYears, days=nDays, CRWCorrelation=CRWCorrelation, nSurveyRoutes=nSurveyRoutes)$newInstance(context=context, isTest=T)
    else if (scenario == "F") stop("unsupported")
    else stop("unsupported")
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
    
    if (plot) {
      boundaryDF <- study$studyArea$toGGDF()
      
      p <- ggplot(boundaryDF, aes(long, lat, group=group)) + geom_path() +
        geom_point(data=mss$initialPopulation$toGGDF(), aes(x, y, group=NA), shape="+", size=8, colour="red", alpha=0.7) + theme_raster()
      plot(p)
      saveFigure(p, filename=paste("SimulatedInitialLocations-", scenario, "-", study$studyArea$region, ".svg", sep=""), bg="transparent")
      
      p <- ggplot(boundaryDF, aes(long, lat, group=group)) + geom_path() +
        geom_path(data=tracks$toGGDF(), aes(long, lat, group=group, colour=id), size=1, alpha=0.7) + theme_raster()
      plot(p)
      saveFigure(p, filename=paste("SimulatedTracks-", scenario, "-", study$studyArea$region, ".svg", sep=""), bg="transparent")
      
      p <- ggplot(boundaryDF, aes(long, lat, group=group)) + geom_path() +
        geom_path(data=surveyRoutes$toGGDF(), aes(long, lat, group=group), size=1, colour="blue", alpha=0.7) + theme_raster()
      plot(p)
      saveFigure(p, filename=paste("SimulatedSurveyRoutes-", scenario, "-", study$studyArea$region, ".svg", sep=""), bg="transparent")
    }
    
    intersections <- SimulatedIntersections$new(study=study, iteration=iteration)
    intersections$findIntersections(tracks, surveyRoutes,  dimension=1)
    
    if (plot) {
      if (F) {
        tracks2 <- tracks$copy()
        tracks2$tracks$herdSize <- 1
        intersections2 <- SimulatedIntersections$new(study=study, iteration=iteration)
        intersections2$findIntersections(tracks2, surveyRoutes,  dimension=1)
        cbind(intersections$intersections@data, intersections2$intersections@data)
      }

      x <- ddply(tracks$tracks, .(id), function(y) data.frame(size=y$herdSize[1]))
      colnames(x) <- c("Agent id", "Herd size")
      p <- ggplot(x, aes(`Agent id`, `Herd size`)) + geom_bar(stat="identity") +
        scale_y_discrete() + theme_minimal()
      plot(p)
      saveFigure(p, filename=paste("SimulatedHerdSize-", scenario, "-", study$studyArea$region, ".svg", sep=""), bg="transparent")
    }
    
    model <- SimulatedSmoothModelSpatioTemporal$new(study=study, iteration=iteration)
    model$setup(intersections=intersections, params=modelParams)
    #plot(model$getUnscaledMesh()); plot(study$studyArea$boundary, add=T, border="blue")
    model$estimate()
    model$collectEstimates()
    
    populationDensity <- model$getPopulationDensity(getSD=FALSE)
    populationSize <- populationDensity$mean$integrate(volume=SimulationPopulationSize$new(study=study, iteration=iteration))
    if (mss$hasHabitatWeights()) {
      habitatSelection <- tracks$getHabitatPreferences(habitatWeightsTemplate=habitatWeights, nSamples=30, save=FALSE)
      habitatWeights$setHabitatSelectionWeights(habitatSelection)
      habitatWeightsRaster <- habitatWeights$getWeightsRaster(save=FALSE)
      #habitatWeights; plot(habitatWeightsRaster)
      populationDensity$mean$weight(habitatWeightsRaster)
    }
    populationSizeWeighted <- populationDensity$mean$integrate(volume=SimulationPopulationSize$new(study=study, iteration=iteration))
    
    habitatWeights; populationSize; populationSizeWeighted
    size <- rbind(size, data.frame(obs=sum(model$data$intersections), fitted=sum(model$data$fittedMean * model$getObservedOffset()), sizew=populationSizeWeighted$sizeData$Estimated, size=populationSize$sizeData$Estimated, iteration=iteration))
    
    
    #populationSize <- model$getPopulationSize(tracks=tracks, withHabitatWeights=mss$hasHabitatWeights())
    #x <- populationSize$sizeData$Estimated
    #model$switchToMesh()
    #populationSize <- model$getPopulationSize(tracks=tracks, withHabitatWeights=mss$hasHabitatWeights())
    
    #size <- rbind(size, data.frame(obs=sum(model$data$intersections), fitted=sum(model$data$fittedMean * model$getObservedOffset()), sizeobs=x, sizenode=populationSize$sizeData$Estimated, iteration=iteration))
  }
  
  return(size)
}

if (task_id > 0) {
  library(doMC)
  registerDoMC(cores=detectCores())
  
  size <- if (task_id == 1) simulate("A")
  else if (task_id == 2) simulate("B")
  else if (task_id == 4) simulate("D")
  else if (task_id == 5) simulate("E")
  
  print(colMeans(size))
  print(colSDs(size))
}
