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
    
    getValidationDataFileName = function() {
      return(context$getFileName(dir=context$processedDataDirectory, name="ValidationPopulationSize", response="", region=study$studyArea$region))
    },
    
    saveValidationData = function() {
      stop("Unimplemented method.")
    },
    
    loadValidationData = function() {
      load(getValidationDataFileName())
      sizeData <<- merge(sizeData, validationSizeData, all=TRUE, sort=FALSE)
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
      sizeData <<- merge(sizeData, validation[,c("Year",study$response)], all=T)
      colnames(sizeData) <<- c("Year","Estimated","Observed")
      return(invisible(.self))
    }
  )
)
