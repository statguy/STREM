Study <- setRefClass(
  Class = "Study",
  fields = list(
    context = "Context",
    response = "character",
    studyArea = "StudyArea"
  ),
  methods = list(
    getTemplateRaster = function() {#height=600, ext=extent(studyArea$habitat)) {
      #library(sp)
      #library(raster)
      
      #dimXY <- dim(studyArea$habitat)[1:2]
      #aspectRatio <- dimXY[2] / dimXY[1]
      #width <- height * aspectRatio
      #templateRaster <- raster(ext, nrows=height, ncols=width, crs=studyArea$proj4string)
      
      height <- dim(studyArea$habitat)[1] / 100 # Determine scaling automatically
      width <- dim(studyArea$habitat)[2] / 100
      templateRaster <- raster(extent(studyArea$habitat), nrows=height, ncols=width, crs=studyArea$proj4string)
      
      return(templateRaster)
    }
  )
)

SimulationStudy <- setRefClass(
  Class = "SimulationStudy",
  contains = "Study",
  methods = list(
    initialize = function(...) {
      callSuper(...)
      return(invisible(.self))
    },
    
    newInstance = function(context, isTest=F) {
      context <<- context
      studyArea <<- if (isTest) TestStudyArea$new(context=context)$newInstance()
      else FinlandStudyArea$new(context=context)$newInstance()
      return(invisible(.self))
    },
    
    loadTracks = function(iteration) {
      tracks <- SimulatedTracks$new(study=.self, iteration=iteration)$loadTracks()
      return(tracks)
    },
    
    loadIntersections = function(iteration) {
      intersections <- SimulatedIntersections$new(study=.self, iteration=iteration)$loadIntersections()
      return(intersections)
    },
    
    loadEstimates = function(iteration) {
      estimates <- SimulatedSmoothModel$new(study=.self, iteration=iteration)$loadEstimates()
      return(estimates)
    },
    
    loadPopulationSize = function(iteration) {
      populationSize <- SimulationPopulationSize$new(study=.self, iteration=iteration)$loadPopulationSize()
      return(populationSize)
    },
    
    #loadTracksCollection = function() {
    #  tracks <- SimulatedTracksCollection$new(study=.self)
    #  tracks$loadTracks()
    #  return(tracks)
    #},
    
    #loadIntersectionsCollection = function() {
    #  intersections <- SimulatedIntersectionsCollection$new(study=.self)
    #  intersections$loadIntersections()
    #  return(intersections)
    #},
    
    loadSurveyRoutes = function(n=800, random=TRUE) {
      surveyRoutes <- if (random) FinlandRandomWTCSurveyRoutes$new(study=.self)$newInstance(n=800)
      else FinlandWTCSurveyRoutes(study=.self)$newInstance()
      return(surveyRoutes)
    }    
  )
)

FinlandWTCStudy <- setRefClass(
  Class = "FinlandWTCStudy",
  contains = "Study",
  fields = list(
  ),
  methods = list(
    initialize = function(context, ...) {
      callSuper(context=context, ...)
      studyArea <<- FinlandStudyArea$new(context=context)$newInstance()
      return(invisible(.self))
    },
    
    preprocessResponse = function(response, maxDuration=1, cacheCovariates=TRUE, findHabitatWeights=TRUE, fmiApiKey) {
      library(CNPCluster)
      
      cnpClusterStartLocal()
      
      response <<- response
      
      intersections <- FinlandWTCIntersections$new(study=.self, maxDuration=maxDuration)
      tracks <- FinlandWTCTracks$new(study=.self)
      habitatWeights <- CORINEHabitatWeights$new(study=.self)

      intersections$saveIntersections()
      intersections$saveCovariates(intersections$intersections, cache=cacheCovariates, fmiApiKey=fmiApiKey)
      
      tracks$saveTracks()
      
      if (findHabitatWeights) {
        habitatSelection <- tracks$getHabitatPreferences(habitatWeightsTemplate=habitatWeights, nSamples=30, save=T)
        habitatWeights <- CORINEHabitatWeights$new(study=.self)$setHabitatSelectionWeights(habitatSelection)
        habitatWeights$getWeightsRaster(save=TRUE)
      }
      
      cnpClusterStopLocal()
    },
    
    preprocess = function(cacheCovariates=TRUE, findHabitatWeights=TRUE, fmiApiKey) {
      preprocessResponse(response="canis.lupus", findHabitatWeights=findHabitatWeights, cacheCovariates=cacheCovariates, fmiApiKey=fmiApiKey)
      preprocessResponse(response="lynx.lynx", findHabitatWeights=findHabitatWeights, cacheCovariates=FALSE, fmiApiKey=fmiApiKey)
      preprocessResponse(response="rangifer.tarandus.fennicus", findHabitatWeights=findHabitatWeights, cacheCovariates=FALSE, fmiApiKey=fmiApiKey)
      return(invisible(.self))
    },
    
    loadTracks = function() {
      return(FinlandWTCTracks$new(study=.self)$loadTracks())
    },
    
    loadSurveyRoutes = function() {
      return(FinlandWTCSurveyRoutes$new(study=.self)$newInstance())
    },
    
    loadIntersections = function() {
      intersections <- FinlandWTCIntersections$new(study=.self)$loadIntersections()
      intersections$loadCovariates()
      return(intersections)
    },
    
    loadHabitatSelection = function() {
      habitatSelection <- HabitatSelection$new(study=.self)$loadHabitatSelection()
      return(habitatSelection)
    },
    
    loadHabitatWeights = function() {
      habitatSelection <- loadHabitatSelection()
      habitatWeights <- CORINEHabitatWeights$new(study=.self)$setHabitatSelectionWeights(habitatSelection=habitatSelection)
      return(habitatWeights)
    },
    
    loadHabitatWeightsRaster = function() {
      return(CORINEHabitatWeights$new(study=.self)$loadWeightsRaster())
    },
    
    estimate = function(test=FALSE, quick=FALSE) {
      meshParams <- if (quick) list(maxEdge=c(.2e6, .4e6), cutOff=.1e6, coordsScale=1e-6)
      else list(maxEdge=c(.05e6, .15e6), cutOff=.02e6, coordsScale=1e-6)
      
      intersections <- loadIntersections()
      intersections$intersections$distance <- 1
            
      model <- FinlandSmoothModel$new(study=.self)
      model$setup(intersections=intersections, meshParams=meshParams)
      xyt <- model$associateMeshLocationsWithDate()
      #model$saveCovariates(xyt, impute=FALSE, save=FALSE)
      model$saveCovariates(xyt, impute=TRUE, save=TRUE)
      
      if (!test) model$estimate(save=T, fileName=model$getEstimatesFileName())
      
      return(model)
    },
        
    loadEstimates = function() {
      estimates <- FinlandSmoothModel$new(study=.self)$loadEstimates()
      estimates$loadCovariates()   
      return(estimates)
    },

    #predictDistances = function(predictCovariates) {
    #  intervals <- FinlandMovementSampleIntervals$new(study=.self)
    #  distances <- intervals$predictDistances(predictCovariates=predictCovariates)
    #  return(distances)
    #},
    
    getSampleIntervals = function() {
      intervals <- FinlandMovementSampleIntervals$new(study=.self)
      return(intervals)
    },
    
    collectEstimates = function(withDistanceWeights=TRUE) {
      intervals <- getSampleIntervals()
      
      message("Predictions for survey routes...")
      intersections <- loadIntersections()
      distancesAtSurveyRoutes <- if (withDistanceWeights)
        subset(intervals$predictDistances(predictCovariates=intersections$covariates), Variable=="Predicted", select="Value", drop=TRUE)
      else
        mean(loadTracks()$getDistances(), na.rm=T)
      
      message("Predictions for mesh nodes...")
      estimates <- loadEstimates()
      distancesAtNodes <- if (withDistanceWeights)
        subset(intervals$predictDistances(predictCovariates=estimates$covariates), Variable=="Predicted", select="Value", drop=TRUE)
      else
        mean(loadTracks()$getDistances(), na.rm=T)
      
      estimates$collectEstimates(weightsAtSurveyRoutes=distancesAtSurveyRoutes, weightsAtNodes=distancesAtNodes, quick=TRUE)
      
      return(estimates)
    },
    
    getPopulationDensity = function(withHabitatWeights=TRUE, withDistanceWeights=TRUE, saveDensityPlots=FALSE) {
      estimates <- collectEstimates(withDistanceWeights=withDistanceWeights)
      
      habitatWeights <- if (withHabitatWeights) loadHabitatWeightsRaster() else HabitatWeights$new(study=study)$getWeightsRaster()
      populationDensity <- estimates$getPopulationDensity(templateRaster=habitatWeights, getSD=FALSE)
      populationDensity$mean$weight(habitatWeights)
      
      if (saveDensityPlots) {
        populationDensity$mean$animate(name=estimates$modelName)
        # TODO: SD
      }
      
      return(populationDensity)
    },
    
    getPopulationSize = function(withHabitatWeights=TRUE, withDistanceWeights=TRUE, saveDensityPlots=FALSE) {
      populationDensity <- getPopulationDensity(withHabitatWeights=withHabitatWeights, withDistanceWeights=withDistanceWeights, saveDensityPlots=saveDensityPlots)
      populationSize <- populationDensity$mean$integrate(volume=FinlandPopulationSize$new(study=study))
      return(populationSize)
    },
    
    show = function() {
      cat("Response: ", response, "\n")
      cat("Study region: ", studyArea$region, "\n")
      loadHabitatSelection()$show()
      loadHabitatWeights()$show()
      return(invisible(.self))
    }
  )
)
