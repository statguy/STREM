HabitatSelection <- setRefClass(
  Class = "HabitatSelection",
  fields = list(
    study = "Study",
    nullModelUsage = "data.frame",
    realizedUsage = "data.frame",
    relativeUsage = "data.frame",
    nCollaredIndividuals = "integer",
    maxTracks = "integer"
  ),
  methods = list(
    initialize = function(nCollaredIndividuals=as.integer(30), maxTracks=as.integer(10000), ...) {
      callSuper(...)
      nCollaredIndividuals <<- nCollaredIndividuals
      maxTracks <<- maxTracks
    },
    
    randomizeSteps = function(movements, location, nSamples) {
      randomizedSteps <- movements[sample(1:nrow(movements), min(nSamples, nrow(movements))), c("dx","dy")]
      randomizedLocations <- cbind(x=location$x + randomizedSteps$dx, y=location$y + randomizedSteps$dy)
      return(randomizedLocations)
    },
    
    getNullModelMovementHabitatDistributions = function(movements, habitatWeightsTemplate, nSamples) {
      library(plyr)
      library(raster)
      library(sp)
      
      movements <- movements[!(is.na(movements$x) | is.na(movements$y) | is.na(movements$dx) | is.na(movements$dy)),]      
      
      p <- dlply(movements, .(id, burst), function(m, habitat, habitatWeightsTemplate, nSamples) {
        message("Processing burst = ", m$burst[1], ", n = ", nrow(m), " for potential movements...")
        if (nrow(m) < 50) {
          message("Aborted, must have at least 50 movement vectors.")
          return(NULL)
        }
        
        locations <- m[,c("x","y")]
        p <- habitatWeightsTemplate$getHabitatFrequencies(c())
        
        for (i in 1:nrow(locations)) {
          location <- locations[i,]
          randomizedLocations <- randomizeSteps(movements=m, location=location, nSamples=nSamples)
          habitatSample <- raster::extract(habitat, SpatialPoints(randomizedLocations, CRS(projection(habitat))))    
          p <- p + habitatWeightsTemplate$getHabitatFrequencies(habitatSample)
        }
        
        return(prop.table(p))
      }, habitat=study$studyArea$habitat, habitatWeightsTemplate=habitatWeightsTemplate, nSamples=nSamples, .parallel=TRUE)
      
      x <- do.call(rbind, p)
      return(as.data.frame(x))
    },
    
    getRealizedMovementHabitatDistributions = function(movements, habitatWeightsTemplate) {
      library(dplyr)
      library(plyr)
      library(raster)
      
      movements <- movements[!(is.na(movements$x) | is.na(movements$y) | is.na(movements$dx) | is.na(movements$dy)),]
      if (nrow(movements) > maxTracks) movements <- sample_n(movements, maxTracks)
      
      p <- dlply(movements, .(id, burst), function(m, habitat, habitatWeightsTemplate) {
        message("Processing burst = ", m$burst[1], " n = ", nrow(m), " for actual movements...")
        if (nrow(m) < 50) {
          message("Aborted, must have at least 50 movement vectors.")
          return(NULL)
        }
                
        locations <- m[,c("x","y")]
        habitatSample <- raster::extract(habitat, SpatialPoints(locations, CRS(projection(habitat))))  
        p <- prop.table(habitatWeightsTemplate$getHabitatFrequencies(habitatSample))
        return(p)
      }, habitat=study$studyArea$habitat, habitatWeightsTemplate=habitatWeightsTemplate, .parallel=TRUE, .paropts=list(.packages=c("raster","sp")))
      
      x <- do.call(rbind, p)
      return(as.data.frame(x))
    },

    getHabitatPreferences = function(intervals, habitatWeightsTemplate, nSamples=10, save=FALSE) {
      if (!inherits(intervals, "MovementSampleIntervals"))
        stop("Argument 'intervals' must be of class 'MovementSampleIntervals'.")
      tracksDF <- intervals$getSampleIntervals()
            
      message("Estimating potential habitat usage...")
      nullModelUsage <<- getNullModelMovementHabitatDistributions(movements=tracksDF, habitatWeightsTemplate=habitatWeightsTemplate, nSamples=nSamples)
      message("Counting actual habitat usage...")
      realizedUsage <<- getRealizedMovementHabitatDistributions(movements=tracksDF, habitatWeightsTemplate=habitatWeightsTemplate)
      relativeUsage <<- as.data.frame(as.list(colSums(realizedUsage) / colSums(nullModelUsage)))
      #relativeUsage95 <<- 
      
      if (save) saveHabitatSelection()
      
      return(invisible(.self))
    },
    
    getHabitatSelectionFileName = function() {
      return(study$context$getFileName(dir=study$context$resultDataDirectory, name="HabitatWeights", response=study$response, region=study$studyArea$region))
    },
    
    saveHabitatSelection = function() {
      fileName <- getHabitatSelectionFileName()
      message("Saving habitat preferences to ", fileName, "...")
      save(nullModelUsage, realizedUsage, relativeUsage, file=fileName)
    },
    
    loadHabitatSelection = function() {
      load(file=getHabitatSelectionFileName(), envir=as.environment(.self))
      return(invisible(.self))
    },
    
    show = function() {
      cat("Habitat null model usage:\n")
      print(colMeans(nullModelUsage))
      cat("Habitat realized usage:\n")
      print(colMeans(realizedUsage))
      cat("Habitat relative usage:\n")
      print(colMeans(relativeUsage))
    }
  )
)

SimulationHabitatSelection <- setRefClass(
  "SimulationHabitatSelection",
  contains = "HabitatSelection",
  fields = list(
    iteration = "integer"  
  ),
  methods = list(
    getHabitatSelectionFileName = function() {
      if (inherits(study, "undefinedField") | length(iteration) == 0)
        stop("Provide study and iteration parameters.")
      return(study$context$getLongFileName(study$context$scratchDirectory, name="HabitatWeights", response=study$response, region=study$studyArea$region, tag=iteration))
    }
  )
)

WTCHabitatSelection <- setRefClass(
  "WTCHabitatSelection",
  contains = "HabitatSelection",
  methods = list(
  )
)
