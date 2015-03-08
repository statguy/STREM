PopulationSize <- setRefClass(
  Class = "PopulationSize",
  fields = list(
    study = "Study",
    sizeData = "data.frame",
    modelName = "character"
  ),
  methods = list(
    initialize = function(...) {
      callSuper(sizeData=data.frame(), ...)
      return(invisible(.self))
    },
    
    addYearSize = function(year, size) {
      sizeData <<- rbind(sizeData, data.frame(Year=year, Estimated=size))
      return(invisible(.self))
    },
    
    setValidationSizeData = function(validationSizeData) {
      sizeData <<- if (nrow(sizeData) == 0) validationSizeData
      else merge(sizeData, validationSizeData, all=TRUE, sort=FALSE)
    },
    
    getValidationDataFileName = function() {
      return(study$context$getFileName(dir=study$context$processedDataDirectory, name="ValidationPopulationSize", response="", region=study$studyArea$region))
    },
    
    saveValidationData = function() {
      stop("Unimplemented method.")
    },
    
    loadValidationData = function(tracks, fileName=getValidationDataFileName()) {
      load(fileName)
      setValidationSizeData(validationSizeData)
      return(invisible(.self))
    },
        
    plotPopulationSize = function() {
      library(ggplot2)
      library(reshape2)
      
      melted <- melt(sizeData, id.vars="Year")
      p <- ggplot(melted, aes(x=Year, y=value, colour=variable, group=variable)) + geom_line() +
        xlab("Year") + ylab("Population size")
      print(p)
      
      return(invisible(.self))
    },
    
    getPopulationSizeFileName = function() {
      if (length(modelName) == 0) stop("Provide modelName parameter.")
      return(study$context$getLongFileName(dir=study$context$resultDataDirectory, name="PopulationSize", response=study$response, region=study$studyArea$region, tag=modelName))
    },
    
    savePopulationSize = function(fileName=getPopulationSizeFileName()) {
      save(sizeData, iteration, file=fileName)
      return(invisible(.self))
    },
    
    loadPopulationSize = function(fileName=getPopulationSizeFileName()) {
      load(fileName, envir=as.environment(.self))
      sizeData$Year <<- as.integer(as.character(sizeData$Year))
      return(invisible(.self))
    },
    
    match = function(model=Validation ~ -1 + Estimated) {
      result <- lm(model, data=sizeData)
      return(result)
    },
    
    show = function() {
      print(sizeData)
      return(invisible(.self))
    },
    
    getPopulationSize = function(density, year, location, habitatWeights, loadValidationData=TRUE) {
      x <- data.frame(density=density, year=year)
      if (missing(location)) {
        x <- ddply(x, .(year), function(x) data.frame(density=mean(x$density), year=x$year[1]))
      }
      x$size <- x$density * study$studyArea$boundary@polygons[[1]]@area
      sizeData <<- data.frame(Year=x$year, Estimated=x$size)
      if (loadValidationData) .self$loadValidationData()
      return(invisible(.self))
    }
    
    
    #getPopulationSize = function(populationDensity, habitatWeights=1, loadHabitatWeights=TRUE, loadValidationData=TRUE) {
    #  if (loadHabitatWeights) habitatWeights <- study$loadHabitatWeightsRaster()
    #  populationDensity$integrate(volume=.self, weights=habitatWeights)      
    #  if (loadValidationData) .self$loadValidationData()
    #  return(invisible(.self))
    #}
  )
)

FinlandPopulationSize <- setRefClass(
  Class = "FinlandPopulationSize",
  contains = "PopulationSize",
  methods = list(    
    saveValidationData = function() {
      # TODO
      return(invisible(.self))
    },
    
    loadValidationData = function(fileName=getValidationDataFileName()) {
      load(fileName)
      validationData <- validation[,c("Year", study$response)]
      colnames(validationData) <- c("Year","Validation")
      setValidationSizeData(validationData)
      return(invisible(.self))
    }
  )
)

SimulationPopulationSize <- setRefClass(
  Class = "SimulationPopulationSize",
  contains = "PopulationSize",
  fields = list(
    iteration = "integer"
  ),
  methods = list(
    getPopulationSizeFileIterations = function() {
      if (inherits(study, "undefinedField") | length(modelName) == 0)
        stop("Provide study and modelName parameters.")
      return(study$context$getIterationIds(dir=study$context$resultDataDirectory, name="PopulationSize", response=study$response, region=study$studyArea$region, tag=paste("(\\d+)", modelName, sep="-")))
    },
    
    getPopulationSizeFileName = function() {
      if (inherits(study, "undefinedField") | length(modelName) == 0 | length(iteration) == 0)
        stop("Provide study, modelName and iteration parameters.")
      return(study$context$getLongFileName(dir=study$context$resultDataDirectory, name="PopulationSize", response=study$response, region=study$studyArea$region, tag=paste(iteration, modelName, sep="-")))
    },
    
    loadValidationData = function(tracks, fileName) {
      if (missing(tracks)) tracks <- study$loadTracks(iteration=iteration, addColumns=FALSE)
      truePopulationSize <- tracks$getTruePopulationSize()
      setValidationSizeData(truePopulationSize)
      return(invisible(.self))
    }
  )
)
