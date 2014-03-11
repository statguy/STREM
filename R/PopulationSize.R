PopulationSize <- setRefClass(
  Class = "PopulationSize",
  fields = list(
    study = "Study",
    sizeData = "data.frame"
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
      sizeData <<- merge(sizeData, validationSizeData, all=TRUE, sort=FALSE)
    },
    
    getValidationDataFileName = function() {
      return(context$getFileName(dir=study$context$processedDataDirectory, name="ValidationPopulationSize", response="", region=study$studyArea$region))
    },
    
    saveValidationData = function() {
      stop("Unimplemented method.")
    },
    
    loadValidationData = function() {
      load(getValidationDataFileName())
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
      return(context$getFileName(dir=study$context$resultDataDirectory, name="PopulationSize", response=study$response, region=study$studyArea$region))
    },
    
    savePopulationSize = function(fileName=getPopulationSizeFileName()) {
      save(sizeData, iteration, file=fileName)
      return(invisible(.self))
    },
    
    loadPopulationSize = function(fileName=getPopulationSizeFileName()) {
      load(fileName, envir=as.environment(.self))
      return(invisible(.self))
    },
    
    show = function() {
      print(sizeData)
      return(invisible(.self))
    }
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
    
    loadValidationData = function() {
      load(getValidationDataFileName())
      setValidationSizeData(validation[,c("Year", study$response)])
      #sizeData <<- merge(sizeData, validation[,c("Year",study$response)], all=T)
      colnames(sizeData) <<- c("Year","Estimated","Observed")
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
    getPopulationSizeFileName = function() {
      return(context$getLongFileName(dir=context$processedDataDirectory, name="PopulationSize", response=study$response, region=study$studyArea$region, tag=iteration))
    },
    
    loadValidationData = function() {
      tracks <- study$loadTracks(iteration=iteration)
      truePopulationSize <- tracks$getTruePopulationSize()
      setValidationSizeData(truePopulationSize)
      return(invisible(.self))
    }
  )
)