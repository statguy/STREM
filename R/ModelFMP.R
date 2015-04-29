FMPModel <- setRefClass(
  Class = "FMPModel",
  #contains = "AggregatedModel",
  contains = "Model",
  fields = list(
  ),
  methods = list(
    setup = function(intersections, params, tag) {
      coordsScale <<- 1
      offsetScale <<- 1000^2
      modelName <<- if (missing(tag)) "FMPModel" else paste("FMPModel", tag, sep="-")
      
      data <<- intersections$getData()
      data$response <<- data$intersections
      locations <<- intersections$getCoordinates() * coordsScale
      years <<- as.integer(sort(unique(data$year)))
      
      #.self$aggregate()
      
      return(invisible(.self))      
    },
    
    saveEstimates = function(fileName=getEstimatesFileName()) {
      message("Saving result to ", fileName, "...")
      save(locations, data, coordsScale, years, offsetScale, file=fileName)
    },
    
    estimate = function(save=FALSE, fileName=getEstimatesFileName(), verbose=TRUE) {
      if (save) saveEstimates(fileName=fileName)
    },
    
    collectEstimates = function(observationWeights=1, predictionWeights=1) {
      index <- 1:nrow(data)
      
      data$fittedMean <<- data$intersections / getObservedOffset() * observationWeights
      data$fittedSD <<- NA
      observedOffset <- getObservedOffset()
      
      message("Fitted values sums all years:")
      message("observed = ", sum(data$intersections))
      message("estimated = ", sum(data$fittedMean * observedOffset))
      
      return(invisible(.self))
    }
  )
)

SimulatedFMPModel <- setRefClass(
  Class = "SimulatedFMPModel",
  contains = "FMPModel",
  fields = list(
    iteration = "integer"
  ),
  methods = list(
    getEstimatesFileIterations = function() {
      if (inherits(study, "undefinedField") | length(modelName) == 0)
        stop("Provide study and modelName parameters.")
      return(study$context$getIterationIds(dir=study$context$scratchDirectory, name=modelName, response=study$response, region=study$studyArea$region, tag="(\\d+)"))
    },
    
    getEstimatesFileName = function() {
      if (inherits(study, "undefinedField") | length(modelName) == 0 | length(iteration) == 0)
        stop("Provide study, modelName and iteration parameters.")
      return(study$context$getLongFileName(study$context$scratchDirectory, name=modelName, response=study$response, region=study$studyArea$region, tag=iteration))
    },
    
    samplePosterior = function(n, index) {
      stop("Unsupported.")
    },
    
    getPopulationSize = function(populationDensity, habitatWeightsRaster=NULL) {
      if (missing(populationDensity))
        stop("Required argument 'populationDensity' missing.")
      if (!inherits(populationDensity, "SpatioTemporalRaster"))
        stop("Argument 'populationDensity' must be of type 'SpatioTemporalRaster'")
      if (!is.null(habitatWeightsRaster)) populationDensity$weight(habitatWeightsRaster)
      
      populationSize <- populationDensity$integrate(volume=SimulationPopulationSize(study=study, iteration=iteration, modelName=modelName))
      populationSize$loadValidationData()
      
      return(invisible(populationSize))
    }
  )
)

FinlandFMPModel <- setRefClass(
  Class = "FinlandFMPModel",
  contains = c("FMPModel", "FinlandCovariates"),
  fields = list(
  ),
  methods = list(
    initialize = function(...) {
      callSuper(covariatesName="FinlandSmoothModelCovariates", ...)
      return(invisible(.self))
    },
    
    predictDistances = function(formula=study$getDistanceCovariatesModel(), intervalH=study$getTrackSampleInterval()) {
      distances <- study$predictDistances(formula=formula, data=covariates, intervalH=intervalH)
      return(distances)
    }
  )
)
