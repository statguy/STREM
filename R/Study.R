Study <- setRefClass(
  Class = "Study",
  fields = list(
    context = "Context",
    response = "character",
    studyArea = "StudyArea"
  ),
  methods = list(
    getTemplateRaster = function() {
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
  fields = list(
    surveyRoutes = "ANY",
    withHabitatWeights = "logical",
    grassLocalTempDir = "character"
  ),
  methods = list(
    initialize = function(...) {
      callSuper(...)
      return(invisible(.self))
    },
    
    setup = function(context, surveyRoutes, withHabitatWeights=F, isTest=F) {
      context <<- context
      withHabitatWeights <<- withHabitatWeights
      if (!missing(surveyRoutes)) surveyRoutes <<- surveyRoutes # WARNING: this can give us recursion for study object references! TODO: better design
      studyArea <<- if (isTest) TestStudyArea$new(context=context)$setup()
      else FinlandStudyArea$new(context=context)$setup()
      return(invisible(.self))
    },
    
    loadTracks = function(iteration, addColumns=TRUE, findTruePopulationSize=TRUE) {
      tracks <- SimulatedTracks$new(study=.self, iteration=iteration)$loadTracks(addColumns=addColumns, findTruePopulationSize=findTruePopulationSize)
      return(tracks)
    },
    
    loadIntersections = function(iteration, predictDistances=TRUE) {
      intersections <- SimulatedIntersections(study=.self, iteration=iteration)$loadIntersections()
      if (predictDistances) {
        tracks <- loadTracks(iteration=iteration)
        intersections$intersections$distance <- tracks$getMeanDistance()
      }
      return(intersections)
    },
    
    getModel = function(modelName, iteration) {
      estimates <- if (startsWith(modelName, "SmoothModel-nbinomial-matern-ar1")) SimulatedSmoothModelSpatioTemporal(study=.self, iteration=iteration)
      else if (startsWith(modelName, "SmoothModel-nbinomial-ar1")) SimulatedSmoothModelTemporal(study=.self, iteration=iteration)
      else if (startsWith(modelName, "SmoothModelMean-nbinomial-ar1")) SimulatedSmoothModelMeanTemporal(study=.self, iteration=iteration)
      else if (startsWith(modelName, "FMPModel")) SimulatedFMPModel(study=.self, iteration=iteration)
      else stop("Invalid model.")
      estimates$modelName <- modelName
      return(estimates)
    },
    
    loadEstimates = function(estimates) {
      estimates$loadEstimates()
      return(estimates)
    },
    
    loadPopulationSize = function(iteration, modelName) {
      populationSize <- SimulationPopulationSize(study=.self, iteration=iteration, modelName=modelName)$loadPopulationSize()$loadValidationData()
      return(populationSize)
    },

    countIntersections = function(iteration, days=1, save=TRUE) {
      tracks <- .self$loadTracks(iteration=iteration)
      intersections <- tracks$countIntersections(surveyRoutes=surveyRoutes, days=days, save=save)
      return(invisible(intersections))
    },
    
    estimate = function(model, params, tag=NULL, predictDistances=FALSE, save=TRUE) {
      if (missing(model)) stop("Missing model argument.")
      if (missing(params)) stop("Missing params argument.")
      
      intersections <- loadIntersections(iteration=model$iteration, predictDistances=predictDistances)
      if (is.null(tag)) model$setup(intersections=intersections, params=params)
      else model$setup(intersections=intersections, params=params, tag=tag)
      model$estimate()
      if (save) model$saveEstimates()
      return(invisible(model))
    },
    
    getHabitatWeights = function(tracks, iteration, save=TRUE, readHabitatIntoMemory=TRUE, asRaster=FALSE) {
      if (withHabitatWeights == FALSE) return(NULL)
      
      habitatWeights <- CORINEHabitatWeights$new(study=.self, iteration=iteration)
      habitatPreferences <- HabitatSelection$new(study=.self, iteration=iteration)
      fileName <- habitatPreferences$getHabitatSelectionFileName()
      habitatSelection <- if (file.exists(fileName)) habitatPreferences$loadHabitatSelection()
      else {      
        message("Habitat preferences file not found: ", fileName)
        message("Estimating...")
        if (missing(tracks)) tracks <- loadTracks(iteration=iteration)
        if (readHabitatIntoMemory)
          study$studyArea$readRasterIntoMemory()
        tracks$getHabitatPreferences(habitatWeightsTemplate=habitatWeights, nSamples=30, save=save)
      }
            
      habitatWeights$setHabitatSelectionWeights(habitatSelection)
      if (asRaster) {
        r <- if (length(grassLocalTempDir) != 0) habitatWeights$getWeightsRaster(save=save, grassLocalTempDir=grassLocalTempDir, iteration=iteration)
        else habitatWeights$getWeightsRaster(save=save)
        return(r)
      }
      return(habitatWeights)
    },
        
    getPopulationSize = function(estimates, index, readHabitatIntoMemory=TRUE, loadValidationData=TRUE, save=TRUE, .parallel=TRUE) {
      if (withHabitatWeights) {
        habitatWeights <- getHabitatWeights(iteration=estimates$iteration, readHabitatIntoMemory=readHabitatIntoMemory)
        populationDensity <- estimates$getPopulationDensity(habitatWeights=habitatWeights, index=index, .parallel=.parallel)
        habitatWeightsRaster <- if (length(grassLocalTempDir) != 0)
          habitatWeights$getWeightsRaster(save=save, grassLocalTempDir=grassLocalTempDir)
        else habitatWeights$getWeightsRaster(save=save)
        populationSize <- estimates$getPopulationSize(populationDensity, habitatWeightsRaster=habitatWeightsRaster)
      }
      else {
        populationSize <- SimulationPopulationSize$new(study=.self, modelName=modelName, iteration=estimates$iteration)
        x <- estimates$getDensityEstimates(index=index)
        populationSize$getPopulationSize(x$density, x$year, loadValidationData=loadValidationData)
      }
      
      if (save) populationSize$savePopulationSize()
      
      return(invisible(populationSize))
    },
    
    getPopulationSize2 = function(modelName, iteration, readHabitatIntoMemory=TRUE, save=TRUE) {
      estimates <- getModel(modelName=modelName, iteration=iteration)
      estimates$loadEstimates()
      estimates$collectEstimates()
      populationSize <- getPopulationSize(estimates, readHabitatIntoMemory=readHabitatIntoMemory, save=save)
      return(invisible(populationSize))
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
      #studyArea <<- FinlandStudyArea$new(context=context)$setup(tolerance=0.001)
      studyArea <<- FinlandStudyArea$new(context=context)$setup()
      return(invisible(.self))
    },
    
    getModel = function(modelName) {
      estimates <- if (modelName == "SmoothModel-nbinomial-matern-ar1") FinlandSmoothModelSpatioTemporal(study=.self)
      else if (modelName == "SmoothModel-nbinomial-ar1") FinlandSmoothModelTemporal(study=.self)
      else if (modelName == "SmoothModel-nbinomial-rw2") FinlandSmoothModelTemporal(study=.self)
      else if (modelName == "SmoothModelMean-nbinomial-ar1") FinlandSmoothModelMeanTemporal(study=.self)
      else if (modelName == "SmoothModelMean-nbinomial-rw2") FinlandSmoothModelMeanTemporal(study=.self)
      else if (modelName == "FMPModel") FinlandFMPModel(study=.self)
      else stop("Invalid model.")
      estimates$modelName <- modelName
      return(estimates)
    },
    
    getModelParams = function(modelName) {
      modelParams <- if (modelName == "SmoothModel-nbinomial-matern-ar1")
        list(family="nbinomial", offsetScale=1000^2, meshParams=.self$studyArea$getMesh(), timeModel="ar1")
      else if (modelName == "SmoothModel-nbinomial-ar1" || modelName == "SmoothModelMean-nbinomial-ar1")
        list(family="nbinomial", offsetScale=1000^2, timeModel="ar1")
      else if (modelName == "SmoothModel-nbinomial-rw2" || modelName == "SmoothModelMean-nbinomial-rw2")
        list(family="nbinomial", offsetScale=1000^2, timeModel="rw2")
      else if (modelName == "FMPModel") NULL
      else stop("Invalid model.")
      return(modelParams)
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
      tracks$saveMetadata()
      
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
      return(FinlandWTCSurveyRoutes$new()$loadSurveyRoutes(context=.self$context, findLengths=findLengths))
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
      habitatSelection <- WTCHabitatSelection$new(study=.self)$loadHabitatSelection()
      return(habitatSelection)
    },
    
    loadHabitatWeights = function() {
      habitatSelection <- loadHabitatSelection()
      habitatWeights <- CORINEHabitatWeights$new(study=.self)$setHabitatSelectionWeights(habitatSelection=habitatSelection)
      return(habitatWeights)
    },
    
    loadHabitatWeightsRaster = function(maskBoundary=TRUE) {
      habitatWeightsRaster <- CORINEHabitatWeights$new(study=.self)$loadWeightsRaster()
      if (maskBoundary)
        habitatWeightsRaster <- mask(habitatWeightsRaster, study$studyArea$boundary)
      return(habitatWeightsRaster)
    },
    
    estimate = function(model, params, predictDistances=TRUE, save=FALSE, test=FALSE) {
      if (missing(model)) stop("Missing model argument.")
      if (missing(params)) stop("Missing params argument.")
      
      intersections <- loadIntersections()
      if (predictDistances) intersections$predictDistances()
      else intersections$intersections$distance <- tracks$getMeanDistance()
      
      model$setup(intersections=intersections, params=params)
      if (!test) model$estimate(save=save, fileName=model$getEstimatesFileName())
      
      return(model)      
    },
    
    loadEstimates = function(estimates) {
      if (missing(estimates)) stop("Missing estimates argument.")
      estimates$loadEstimates()$loadCovariates()   
      return(estimates)
    },

    #getSampleIntervals = function() {
    #  return(FinlandMovementSampleIntervals$new(study=.self))
    #},
    
    #loadSampleIntervals = function() {
    #  return(getSampleIntervals()$loadSampleIntervals())
    #},
    
    #getDistanceCovariatesModel = function() {
    #  return(distanceCovariatesModel);
    #},
    
    #getTrackSampleInterval = function() {
    #  return(trackSampleInterval);
    #},
    
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
    
    collectEstimates = function(model) {
      estimates <- loadEstimates(model)
      estimates$collectEstimates()
      return(estimates)
    },
    
    getPopulationDensity = function(modelName, withHabitatWeights=TRUE, getSD=FALSE, saveDensityPlots) {
      habitatWeightsRaster <- if (withHabitatWeights) loadHabitatWeightsRaster() else 1
      model <- getModel(modelName=modelName)  
      model$offsetScale <- 1000^2 # quickfix, remove when not needed anymore
      model$loadEstimates()
      model$collectEstimates()
      populationDensity <- model$getPopulationDensity(templateRaster=habitatWeightsRaster, getSD=getSD)
      
      #habitatWeights <- if (withHabitatWeights) loadHabitatWeightsRaster() else HabitatWeights$new(study=study)$getWeightsRaster()
      
      if (saveDensityPlots) {
        populationDensity$mean$animate(name="WeightedPopulationDensity-mean")
        if (getSD) populationDensity$sd$animate(name="WeightedPopulationDensity-sd")
      }
      
      return(populationDensity)
    },
        
    getPopulationSize = function(modelName) {
      #populationDensity <- getPopulationDensity(model=model, withHabitatWeights=withHabitatWeights, saveDensityPlots=saveDensityPlots)
      #populationSize <- populationDensity$mean$integrate(volume=FinlandPopulationSize$new(study=study))
      
      habitatWeightsRaster <- loadHabitatWeightsRaster()
      model <- getModel(modelName=modelName)  
      model$offsetScale <- 1000^2 # quickfix, remove when not needed anymore
      model$loadEstimates()
      model$collectEstimates()
      populationDensity <- model$getPopulationDensity(templateRaster=habitatWeightsRaster, getSD=FALSE)
      populationSize <- FinlandPopulationSize(study=.self, modelName=modelName)$getPopulationSize(populationDensity=populationDensity$mean, habitatWeights=habitatWeightsRaster, loadHabitatWeights=FALSE)
      #populationSize$plotPopulationSize()
      #print(populationSize)
      
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
      studyArea <<- RussiaStudyArea$new(context=context)$setup()
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
    },
    
    loadIntersections = function(predictDistances=TRUE) {
      intersections <- RussiaWTCIntersections$new(study=.self)
      intersections$loadIntersections()
      #tracks <- loadTracks()
      #intersections$intersections$distance <- tracks$getMeanDistance()
      return(intersections)
    }
  )
)

FinlandRussiaWTCStudy <- setRefClass(
  Class = "FinlandRussiaWTCStudy",
  contains = "Study",
  methods = list(
    initialize = function(context, ...) {
      callSuper(context=context, ...)
      studyArea <<- FinlandRussiaStudyArea$new(context=context)$setup()
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
    
    loadIntersections = function(predictDistances=TRUE) {
      intersections <- FinlandRussiaWTCIntersections$new(study=.self)
      intersections$loadIntersections()
      if (predictDistances) {
        tracks <- loadTracks()
        intersections$intersections$distance <- tracks$getMeanDistance()
      }
      return(intersections)
    },
    
    estimate = function(model, params, save=FALSE, test=FALSE) {
      if (missing(model)) stop("Missing model argument.")
      if (missing(params)) stop("Missing params argument.")
      
      intersections <- loadIntersections()
      model$setup(intersections=intersections, params=params)      
      if (!test) model$estimate(save=save, fileName=model$getEstimatesFileName())
      
      return(model)
    },
    
    loadEstimates = function(estimates) {
      if (missing(estimates)) stop("Missing estimates argument.")
      estimates$loadEstimates()
      return(estimates)
    },
    
    collectEstimates = function() {
      estimates <- loadEstimates()
      estimates$collectEstimates()
      return(estimates)
    },
    
    getPopulationDensity = function(saveDensityPlots=FALSE, getSD=FALSE) {
      estimates <- collectEstimates()
      
      habitatWeights <- HabitatWeights$new(study=study)$getWeightsRaster()
      populationDensity <- estimates$getPopulationDensity(templateRaster=habitatWeights, getSD=getSD)
      populationDensity$mean$weight(habitatWeights)
      
      if (saveDensityPlots) {
        populationDensity$mean$animate(name=estimates$modelName)
        if (getSD) populationDensity$sd$animate(name=estimates$modelName)
      }
      
      return(populationDensity)
    },
    
    getPopulationSize = function(saveDensityPlots=FALSE) {
      populationDensity <- getPopulationDensity(saveDensityPlots=saveDensityPlots)
      populationSize <- populationDensity$mean$integrate(volume=FinlandRussiaPopulationSize$new(study=study))
      return(populationSize)
    }
  )
)
