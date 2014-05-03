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
    
    preprocess = function() {
      loadSurveyRoutes(save=TRUE)
      return(invisible(.self))
    },
    
    loadTracks = function(iteration, addColumns=TRUE) {
      tracks <- SimulatedTracks$new(study=.self, iteration=iteration)$loadTracks(addColumns=addColumns)
      return(tracks)
    },
    
    loadIntersections = function(iteration) {
      intersections <- SimulatedIntersections$new(study=.self, iteration=iteration)$loadIntersections()
      tracks <- loadTracks(iteration=iteration)
      intersections$intersections$distance <- tracks$getMeanDistance()
      return(intersections)
    },
    
    loadEstimates = function(iteration) {
      estimates <- SimulatedSmoothModel$new(study=.self, iteration=iteration)$loadEstimates()
      return(estimates)
    },
    
    loadPopulationSize = function(iteration) {
      populationSize <- SimulationPopulationSize$new(study=.self, iteration=iteration)$loadPopulationSize()$loadValidationData()
      return(populationSize)
    },
    
    loadSurveyRoutes = function(n=800, random=TRUE, save=FALSE, findLengths=TRUE) {
      surveyRoutes <- if (random) {
        surveyRoutes <- FinlandRandomWTCSurveyRoutes$new(study=.self)
        if (save) surveyRoutes$randomizeSurveyRoutes(nSurveyRoutes=n, save=TRUE)
        return(surveyRoutes$loadSurveyRoutes(findLengths=findLengths))
      }
      else {
        return(FinlandWTCSurveyRoutes$new(study=.self)$loadSurveyRoutes(findLengths=findLengths))
      }
    },
    
    countIntersections = function(iteration) {
      tracks <- study$loadTracks(iteration=iteration)
      intersections <- tracks$countIntersections()
      return(invisible(intersections))
    },
    
    estimate = function(iteration, meshParams, interceptPriorParams, save=TRUE) {
      intersections <- loadIntersections(iteration=iteration)
      model <- SimulatedSmoothModel$new(study=.self, iteration=iteration)
      model$setup(intersections=intersections, meshParams=meshParams)
      if (!missing(interceptPriorParams)) model$setupInterceptPrior(interceptPriorParams)
      model$estimate()
      if (save) model$saveEstimates()
      return(invisible(model))
    },
    
    getPopulationSize = function(iteration, withHabitatWeights=FALSE, save=TRUE) {
      estimates <- loadEstimates(iteration=iteration)
      estimates$collectEstimates()
      populationSize <- estimates$getPopulationSize(withHabitatWeights=withHabitatWeights)
      if (save) populationSize$savePopulationSize()
      return(invisible(populationSize))
    }
  )
)

SimulationStudySubset <- setRefClass(
  Class = "SimulationStudySubset",
  contains = "SimulationStudy",
  fields = list(
    years = "integer"
  ),
  methods = list(
    loadTracks = function(iteration, addColumns=TRUE) {
      tracks <- callSuper(iteration=iteration, addColumns=addColumns)
      tracks$tracks <- subset(tracks$tracks, year %in% years)
      return(tracks)
    },
    
    loadIntersections = function(iteration) {
      intersections <- callSuper(iteration=iteration)
      intersections$intersections <- subset(intersections$intersections, year %in% years)
      return(intersections)
    }
  )
)

FinlandWTCStudy <- setRefClass(
  Class = "FinlandWTCStudy",
  contains = "Study",
  fields = list(
    distanceCovariatesModel = "formula",
    trackSampleInterval = "numeric"
  ),
  methods = list(
    initialize = function(context, ...) {
      callSuper(context=context, ...)
      studyArea <<- FinlandStudyArea$new(context=context)$newInstance(tolerance=0.001)
      return(invisible(.self))
    },
    
    getPrettyResponse = function(response) {
      x <- if (missing(response)) .self$response else response
      y <- switch (x,
        canis.lupus="Canis lupus",
        lynx.lynx="Lynx lynx",
        rangifer.tarandus.fennicus="Rangifer tarandus fennicus")
      return(y)
    },
    
    preprocessResponse = function(response, maxDuration=1, cacheCovariates=TRUE, findHabitatWeights=TRUE, fmiApiKey) {
      response <<- response
      
      intersections <- FinlandWTCIntersections$new(study=.self, maxDuration=maxDuration)
      tracks <- FinlandWTCTracks$new(study=.self)
      habitatWeights <- CORINEHabitatWeights$new(study=.self)

      intersections$saveIntersections()
      intersections$saveCovariates(intersections$intersections, cache=cacheCovariates, fmiApiKey=fmiApiKey)
            
      tracks$saveTracks()
      
      fitDistanceCorrectionModel(save=TRUE)
      
      if (findHabitatWeights) {
        habitatSelection <- tracks$getHabitatPreferences(habitatWeightsTemplate=habitatWeights, nSamples=30, save=T)
        habitatWeights <- CORINEHabitatWeights$new(study=.self)$setHabitatSelectionWeights(habitatSelection)
        habitatWeights$getWeightsRaster(save=TRUE)
      }
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
    
    loadSurveyRoutes = function(findLengths=TRUE) {
      return(FinlandWTCSurveyRoutes$new(study=.self)$loadSurveyRoutes(findLengths=findLengths))
    },
    
    loadIntersections = function(predictDistances=TRUE) {
      intersections <- FinlandWTCIntersections$new(study=.self)
      intersections$loadIntersections()
      intersections$loadCovariates()
      
      if (predictDistances) intersections$predictDistances()
      else {
        tracks <- loadTracks()
        intersections$intersections$distance <- tracks$getMeanDistance()
      }
      
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
    
    estimate = function(meshParams, interceptPriorParams, predictDistances=TRUE, save=FALSE, test=FALSE) {
      intersections <- loadIntersections()
      if (predictDistances) intersections$predictDistances()
      else intersections$intersections$distance <- tracks$getMeanDistance()
      
      model <- FinlandSmoothModel$new(study=.self)
      model$setup(intersections=intersections, meshParams=meshParams)
      if (!missing(interceptPriorParams))
        model$setupInterceptPrior(priorParams=interceptPriorParams)
      
      if (!test) model$estimate(save=save, fileName=model$getEstimatesFileName())
      
      return(model)
    },
        
    loadEstimates = function() {
      estimates <- FinlandSmoothModel$new(study=.self)$loadEstimates()
      estimates$loadCovariates()   
      return(estimates)
    },

    getSampleIntervals = function() {
      return(FinlandMovementSampleIntervals$new(study=.self))
    },
    
    loadSampleIntervals = function() {
      return(getSampleIntervals()$loadSampleIntervals())
    },
    
    getDistanceCovariatesModel = function() {
      return(distanceCovariatesModel);
    },
    
    getTrackSampleInterval = function() {
      return(trackSampleInterval);
    },
    
    fitDistanceCorrectionModel = function(iterations=1000, chains=1, save=TRUE) {
      tracks <- loadTracks()
      
      sampleIntervals <- getSampleIntervals()
      sampleIntervals$getThinnedTracksSampleIntervals(tracks=tracks)
      sampleIntervals$saveTrackCovariates(save=FALSE)
      sampleIntervals$fit(covariatesFormula=getDistanceCovariatesModel(), iterations=iterations, chains=chains)
      if (save) sampleIntervals$saveSampleIntervals()
      
      return(sampleIntervals)
    },
    
    predictDistances = function(formula, data, intervalH=2) {
      if (is.null(formula) | length(intervalH) == 0)
        stop("Argument missing.")
      sampleIntervals <- loadSampleIntervals()
      fixed_model_matrix <- model.matrix(formula, data)
      distances <- sampleIntervals$predict(fixed_model_matrix=fixed_model_matrix, intervalH=intervalH) * 1000
      return(distances)
    },
    
    collectEstimates = function() {
      estimates <- loadEstimates()
      estimates$collectEstimates()
      return(estimates)
    },
    
    getPopulationDensity = function(withHabitatWeights=TRUE, saveDensityPlots=FALSE, getSD=FALSE) {
      estimates <- collectEstimates()
      
      habitatWeights <- if (withHabitatWeights) loadHabitatWeightsRaster() else HabitatWeights$new(study=study)$getWeightsRaster()
      populationDensity <- estimates$getPopulationDensity(templateRaster=habitatWeights, getSD=getSD)
      populationDensity$mean$weight(habitatWeights)
      
      if (saveDensityPlots) {
        populationDensity$mean$animate(name=estimates$modelName)
        if (getSD) populationDensity$sd$animate(name=estimates$modelName)
      }
      
      return(populationDensity)
    },
    
    getPopulationSize = function(withHabitatWeights=TRUE, saveDensityPlots=FALSE) {
      populationDensity <- getPopulationDensity(withHabitatWeights=withHabitatWeights, saveDensityPlots=saveDensityPlots)
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

RussiaWTCStudy <- setRefClass(
  Class = "RussiaWTCStudy",
  contains = "Study",
  methods = list(
    initialize = function(context, ...) {
      callSuper(context=context, ...)
      studyArea <<- RussiaStudyArea$new(context=context)$newInstance()
      return(invisible(.self))
    },
    
    preprocessResponse = function(response) {
      response <<- response
      intersections <- RussiaWTCIntersections$new(study=.self)
      intersections$saveIntersections()
    },
    
    preprocess = function() {
      preprocessResponse(response="canis.lupus")
      preprocessResponse(response="lynx.lynx")
      preprocessResponse(response="rangifer.tarandus.fennicus")
      return(invisible(.self))
    }
  )
)

FinlandRussiaWTCStudy <- setRefClass(
  Class = "FinlandRussiaWTCStudy",
  contains = "Study",
  methods = list(
    initialize = function(context, ...) {
      callSuper(context=context, ...)
      studyArea <<- FinlandRussiaStudyArea$new(context=context)$newInstance()
      return(invisible(.self))
    },
    
    preprocessResponse = function(response) {
      response <<- response
      intersections <- FinlandRussiaWTCIntersections$new(study=.self)
      intersections$saveIntersections()
    },
    
    preprocess = function() {
      preprocessResponse(response="canis.lupus")
      preprocessResponse(response="lynx.lynx")
      preprocessResponse(response="rangifer.tarandus.fennicus")
      return(invisible(.self))
    },
    
    loadTracks = function() {
      finlandStudy <- FinlandWTCStudy$new(context=.self$context, response=response)
      return(FinlandWTCTracks$new(study=finlandStudy)$loadTracks())
    },
    
    loadIntersections = function() {
      intersections <- FinlandRussiaWTCIntersections$new(study=.self)
      intersections$loadIntersections()
      tracks <- loadTracks()
      intersections$intersections$distance <- tracks$getMeanDistance()
      return(intersections)
    },
    
    estimate = function(meshParams, interceptPriorParams, save=FALSE, test=FALSE) {
      intersections <- loadIntersections()
      
      model <- FinlandRussiaSmoothModel$new(study=.self)
      model$setup(intersections=intersections, meshParams=meshParams)
      if (!missing(interceptPriorParams))
        model$setupInterceptPrior(priorParams=interceptPriorParams)
      
      if (!test) model$estimate(save=save, fileName=model$getEstimatesFileName())
      
      return(model)
    }
  )
)
