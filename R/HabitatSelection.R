HabitatSelection <- setRefClass(
  Class = "HabitatSelection",
  fields = list(
    study = "Study",
    tracks = "Tracks",
    expected = "data.frame",
    realized = "data.frame",
    relative = "data.frame"
  ),
  methods = list(
    getExpectedHabitatDistributions = function(movements, habitatWeightsTemplate, nSamples) {
      library(plyr)
      library(raster)
      
      movements <- movements[!(is.na(movements$x) | is.na(movements$y) | is.na(movements$dx) | is.na(movements$dy)),]      
      
      p <- dlply(movements, .(id, burst), function(m, habitat, habitatWeightsTemplate, nSamples) {
        message("Burst = ", m$burst[1], ", n = ", nrow(m))
        
        locations <- m[,c("x","y")]
        p <- habitatWeightsTemplate$getHabitatFrequencies(c())
        
        for (i in 1:nrow(locations)) {
          location <- locations[i,]
          randomizedSteps <- m[sample(1:nrow(m), min(nSamples, nrow(m))), c("dx","dy")]
          randomizedLocations <- cbind(x=location$x + randomizedSteps$dx, y=location$y + randomizedSteps$dy)
          habitatSample <- extract(habitat, SpatialPoints(randomizedLocations, CRS(projection(habitat))))    
          p <- p + habitatWeightsTemplate$getHabitatFrequencies(habitatSample)
        }
        
        return(prop.table(p))
      }, habitat=study$studyArea$habitat, habitatWeightsTemplate=habitatWeightsTemplate, nSamples=nSamples, .parallel=TRUE, .paropts=list(.packages=c("raster","sp")))
      
      x <- do.call(rbind, p)
      return(as.data.frame(x))
    },
    
    getRealizedMovementHabitatDistributions = function(movements, habitatWeightsTemplate) {
      library(plyr)
      library(raster)
      
      movements <- movements[!(is.na(movements$x) | is.na(movements$y) | is.na(movements$dx) | is.na(movements$dy)),]      
      
      p <- dlply(movements, .(id, burst), function(m, habitat, habitatWeightsTemplate) {
        message("Processing burst = ", m$burst[1], " n = ", nrow(m), "...")
        
        locations <- m[,c("x","y")]
        habitatSample <- extract(habitat, SpatialPoints(locations, CRS(projection(habitat))))  
        p <- prop.table(habitatWeightsTemplate$getHabitatFrequencies(habitatSample))
        return(p)
      }, habitat=study$studyArea$habitat, habitatWeightsTemplate=habitatWeightsTemplate, .parallel=TRUE, .paropts=list(.packages=c("raster","sp")))
      
      x <- do.call(rbind, p)
      return(as.data.frame(x))
    },
    
    getHabitatPreferences = function(tracks, habitatWeightsTemplate, nSamples=30, save=FALSE) {
      intervals <- tracks$getSampleIntervals()
      maxIntervalH <- as.numeric(names(which.max(table(intervals$intervals$intervalH))))
      maxIntervals <- subset(intervals$intervals, intervalH == maxIntervalH)
      tracksDF <- ld(tracks$tracks)
      tracksDF <- subset(tracksDF, burst %in% maxIntervals$burst & id %in% maxIntervals$id)
      
      #habitatWeightsTemplate <- CORINEHabitatWeights$new()
      expected <<- getExpectedMovementHabitatDistributions(movements=tracksDF, habitatWeightsTemplate=habitatWeightsTemplate, nSamples=nSamples)
      realized <<- getRealizedMovementHabitatDistributions(movements=tracksDF, habitatWeightsTemplate=habitatWeightsTemplate)
      relative <<- colSums(realizedMovementHabitatUsage) / colSums(expectedMovementHabitatUsage)
      
      #forestRelativeMovementHabitatUsage <- relativeMovementHabitatUsage / relativeMovementHabitatUsage["Forest"]
    }
    
    
  )
)
