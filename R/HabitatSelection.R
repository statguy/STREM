HabitatSelection <- setRefClass(
  Class = "HabitatSelection",
  fields = list(
    study = "Study",
    tracks = "Tracks",
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
      
      movements <- movements[!(is.na(movements$x) | is.na(movements$y) | is.na(movements$dx) | is.na(movements$dy)),]      
      
      p <- dlply(movements, .(id, burst), function(m, habitat, habitatWeightsTemplate, nSamples) {
        library(raster)
        library(sp)
        message("Burst = ", m$burst[1], ", n = ", nrow(m))
        
        locations <- m[,c("x","y")]
        p <- habitatWeightsTemplate$getHabitatFrequencies(c())
        
        for (i in 1:nrow(locations)) {
          location <- locations[i,]
          randomizedLocations <- randomizeSteps(movements=m, location=location, nSamples=nSamples)
          #randomizedSteps <- m[sample(1:nrow(m), min(nSamples, nrow(m))), c("dx","dy")]
          #randomizedLocations <- cbind(x=location$x + randomizedSteps$dx, y=location$y + randomizedSteps$dy)
          habitatSample <- raster::extract(habitat, SpatialPoints(randomizedLocations, CRS(projection(habitat))))    
          p <- p + habitatWeightsTemplate$getHabitatFrequencies(habitatSample)
        }
        
        return(prop.table(p))
      }, habitat=study$studyArea$habitat, habitatWeightsTemplate=habitatWeightsTemplate, nSamples=nSamples, .parallel=TRUE)
      
      x <- do.call(rbind, p)
      return(as.data.frame(x))
    },
    
    getRealizedMovementHabitatDistributions = function(movements, habitatWeightsTemplate) {
      library(plyr)
      library(raster)
      library(dplyr)
      
      movements <- movements[!(is.na(movements$x) | is.na(movements$y) | is.na(movements$dx) | is.na(movements$dy)),]
      if (nrow(movements) > maxTracks) movements <- sample_n(movements, maxTracks)
      
      p <- dlply(movements, .(id, burst), function(m, habitat, habitatWeightsTemplate) {
        message("Processing burst = ", m$burst[1], " n = ", nrow(m), "...")
        
        locations <- m[,c("x","y")]
        habitatSample <- raster::extract(habitat, SpatialPoints(locations, CRS(projection(habitat))))  
        p <- prop.table(habitatWeightsTemplate$getHabitatFrequencies(habitatSample))
        return(p)
      }, habitat=study$studyArea$habitat, habitatWeightsTemplate=habitatWeightsTemplate, .parallel=TRUE, .paropts=list(.packages=c("raster","sp")))
      
      x <- do.call(rbind, p)
      return(as.data.frame(x))
    },

    getMovements = function(tracks) {
      tracks$sample(nCollaredIndividuals)
      intervals <- tracks$getSampleIntervals(collaredIndividuals=collaredIndividuals)
      maxIntervalH <- as.numeric(names(which.max(table(intervals$intervals$intervalH))))
      maxIntervals <- subset(intervals$intervals, intervalH == maxIntervalH)

      tracksDF <- if (inherits(tracks$tracks, "ltraj")) ld(tracks$tracks) else tracks$tracks
      tracksDF <- subset(tracksDF, burst %in% maxIntervals$burst & id %in% maxIntervals$individualId)
      
      #library(dplyr)
      #tracks$tracks %.%
      #  filter(burst %in% maxIntervals$burst, id %in% maxIntervals$individualId) %.%
      #  group_by(burst) %.%
      #  mutate(dx=x[2:n()] - x[1:(n()-1)], dy=y[2:n()] - y[1:(n()-1)])
      
      return(tracksDF)
    },
    
    # TODO: Remove the movements when the individuals are not moving
    getHabitatPreferences = function(tracks, habitatWeightsTemplate, nSamples=30, save=FALSE) {
      message("Finding tracks with highest number of samples and constant frequency...")
      tracksDF <- getMovements(tracks)
      
      message("Estimating potential habitat usage...")
      nullModelUsage <<- getNullModelMovementHabitatDistributions(movements=tracksDF, habitatWeightsTemplate=habitatWeightsTemplate, nSamples=nSamples)
      message("Counting actual habitat usage...")
      realizedUsage <<- getRealizedMovementHabitatDistributions(movements=tracksDF, habitatWeightsTemplate=habitatWeightsTemplate)
      relativeUsage <<- as.data.frame(as.list(colSums(realizedUsage) / colSums(nullModelUsage)))
      
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
    
    plotSampleSteps = function(tracks, plot=TRUE, index=1:5) {
      library(sp)
      library(raster)
      library(rasterVis)
      library(grid)
      
      movements <- getMovements(tracks)
      movements <- movements[!(is.na(movements$x) | is.na(movements$y) | is.na(movements$dx) | is.na(movements$dy)),]
      location <- movements[index[length(index)/2],,drop=FALSE]
      path <- movements[index,c("x","y")]
      randomizedLocations <- randomizeSteps(movements=movements, location=location, nSamples=30)
      randomVectors <- adply(randomizedLocations, 1, function(y, x) data.frame(x=x[1], y=x[2], xend=y[1], yend=y[2]), x=location[,c("x","y"),drop=F])
      
      xy <- SpatialPoints(rbind(path, randomizedLocations), proj4string=study$studyArea$proj4string)
      habitat <- crop(study$studyArea$habitat, extent(xy) * 1.1)
      
      p <- gplot(habitat) + geom_raster(aes(fill=as.factor(value))) +
        scale_fill_manual(values=habitat@legend@colortable) +
        coord_equal() + theme_raster() +
        geom_segment(data=randomVectors, aes(x=x, y=y, xend=xend, yend=yend), colour="darkgrey", size=2, arrow=arrow()) +
        geom_point(data=randomVectors, aes(x=xend, y=yend), colour="darkgrey", size=6, alpha=1) +
        geom_segment(data=path, aes(x=x, y=y, xend=c(tail(x,-1), NA), yend=c(tail(y,-1), NA) ), colour="black", size=2, arrow=arrow()) +
        geom_point(data=path, aes(x=x, y=y), colour="black", size=6, alpha=1)
      if (plot) print(p)
      
      return(p)
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
