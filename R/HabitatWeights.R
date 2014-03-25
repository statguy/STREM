HabitatWeights <- setRefClass(
  Class = "HabitatWeights",
  fields = list(
    study = "Study"
  ),
  methods = list(
    initialize = function(...) {
      callSuper(...)
      return(invisible(.self))
    },
    
    classify = function(habitatValues) {
      return(rep(1, times=length(habitatValues)))
    },
    
    getWeights = function(habitatValues) {
      return(rep(1, times=length(habitatValues)))
    },

    getWeightsRasterFileName = function() {
      study$context$getFileName(dir=study$context$resultDataDirectory, name="HabitatWeightsRaster", response=study$response, region=study$studyArea$region, ext="")
    },

    getWeightsRaster = function(save=FALSE) { # TODO: determine aggregation scale automatically
      if (save) stop("Saving unsupported.")
      library(raster)
      weightsRaster <- study$getTemplateRaster()
      weightsRaster[] <- 1
      return(weightsRaster)
    },

    loadWeightsRaster = function() {
      return(raster(getWeightsRasterFileName()))
    },
    
    getHabitatFrequencies = function(habitatValues) {
      stop("Implement this method in a subclass.")
    },
    
    setHabitatSelectionWeights = function(relativeHabitatUsage) {
      stop("Implement this method in a subclass.")
    },
    
    show = function() {
      cat("Habitat weights:\n1\n")
      return(invisible(.self))
    }
  )
)

CORINEHabitatWeights <- setRefClass(
  Class = "CORINEHabitatWeights",
  contains = "HabitatWeights",
  fields = list(
    weights = "data.frame",
    habitatTypes = "integer"
  ),
  methods = list(
    initialize = function(habitatWeightsList=list(Urban=1, Agriculture=1, Forestland=1, Peatland=1, Water=1), ...) {
      callSuper(...)
      
      lastIndex <- 1+13+4+18+6+3
      
      weights <<- data.frame(
        habitat=as.integer(0:255),
        type=as.integer(c(
          rep(0,1),  # unknown
          rep(1,13), # urban
          rep(2,4),  # agriculture
          rep(3,18), # forestland
          rep(4,6),  # peatland
          rep(5,3),  # water
          rep(0,256-lastIndex)) # undefined
        ),
        weight=c(
          rep(0,1), # unknown
          rep(habitatWeightsList$Urban,13),
          rep(habitatWeightsList$Agriculture,4),
          rep(habitatWeightsList$Forestland,18),
          rep(habitatWeightsList$Peatland,6),
          rep(habitatWeightsList$Water,3),
          rep(0,256-lastIndex) # undefined
        )
      )
      
      habitatTypes <<- unique(weights$type[weights$type != 0])
      names(habitatTypes) <<- c("Urban","Agriculture","Forestland","Peatland","Water")
      
      return(invisible(.self))
    },

    setHabitatSelectionWeights = function(habitatSelection) {
      w <- as.list(habitatSelection$relativeUsage / habitatSelection$relativeUsage[["Forestland"]])
      initialize(habitatWeightsList=w)
      return(invisible(.self))
    },
    
    classify = function(habitatValues, na.value=NA) {
      x <- weights$type[match(habitatValues, weights$habitat)]
      x[is.na(x)] <- na.value
      return(x)
    },

    getWeights = function(habitatValues, na.value=NA) {
      x <- weights$weight[match(habitatValues, weights$habitat)]
      x[is.na(x)] <- na.value
      return(x)
    },
    
    getHabitatFrequencies = function(habitatValues) {
      mappedValues <- classify(as.vector(habitatValues))
      mappedValues <- mappedValues[!(is.na(mappedValues))]
      p <- table(factor(mappedValues,
                        levels=habitatTypes,
                        labels=c("Urban", "Agriculture", "Forestland", "Peatland", "Water")))
      return(p)
    },
    
    getWeightsRaster = function(save=FALSE) { # TODO: determine aggregation scale automatically
      library(raster)
      
      if (save) {
        fileName <- getWeightsRasterFileName()
        
        message("Saving habitat weights raster to ", fileName, ".grd...")
        weightsRaster <- aggregate(study$studyArea$habitat,
                                   aggregationScale=100, # TODO: get this from template raster
                                   filename=fileName, overwrite=TRUE,
                                   fun=function(habitatValue, na.rm)
                                     mean(getWeights(habitatValue), na.rm=na.rm),
                                   na.rm=T)
      }
      else {
        weightsRaster <- aggregate(study$studyArea$habitat,
                                   aggregationScale=100,
                                   fun=function(habitatValue, na.rm)
                                     mean(getWeights(habitatValue), na.rm=na.rm), na.rm=T)
      }
      
      return(weightsRaster)
    },
    
    getHabitatSelectionWeights = function() {
      w <- numeric(length(habitatTypes))
      names(w) <- names(habitatTypes)
      for (type in habitatTypes) w[type] <- weights$weight[weights$type == type][1]
      return(w)
    },
    
    show = function() {
      cat("Habitat weights:\n")
      print(getHabitatSelectionWeights())
      return(invisible(.self))
    }
  )
)
