TracksCollection <- setRefClass(
  Class = "TracksCollection",
  fields = list(
    study = "Study",
    tracksList = "list",
    tracksFiles = "list"
  ),
  methods = list(
    initialize = function(...) {
      callSuper(...)
      return(invisible(.self))
    },
  
    getNumberOfIterations = function() return(length(tracksList)),
    
    addTracks = function(tracks) {
      tracksList[[getNumberOfIterations() + 1]] <<- tracks
      return(invisible(.self))
    },
    
    getTracks = function(index) {
      if (getNumberOfIterations() < index) stop("Invalid index = ", index, ".")
      return(tracksList[[index]])
    },
        
    applyTracks = function(fun, ...) {
      return(lapply(X=tracksList, FUN=fun, ...))
    },
    
    loadTracks = function() {
      stop("Unimplemented method.")
    },    
    
    findIntersections = function(surveyRoutes, dimension) {
      stop("Unimplemented method.")
    },
    
    getDistances = function() {
      applyTracks(function(tracks) return(tracks$getDistances()))
      return(invisible(.self))
    }
  )
)

SimulatedTracksCollection <- setRefClass(
  Class = "SimulatedTracksCollection",
  contains = "TracksCollection",
  methods = list(
    loadTracks = function() {
      # TODO: bad design here if directory or file name pattern changes in Tracks (sub)class...
      tracksFiles <<- study$context$listLongFiles(dir=study$context$scratchDirectory, name="Tracks", response=study$response, tag="\\d+", region=study$studyArea$region)
      if (length(tracksFiles) == 0)
        stop("Cannot find any tracks files in ", study$context$scratchDirectory)
      
      #if (onlyFile == FALSE) {
      #  for (iteration in 1:length(tracksFiles)) { # TODO: better to use iteration number rather than number of files
      #    tracks <- SimulatedTracks$new(study=study, iteration=iteration)
      #    tracks$loadTracks()
      #    addTracks(tracks)
      #  }
      #}
    },
    
    DEPRECATED.randomizeObservationDays = function(minYday=0, maxYday=60) {
      tracksList <<- apply(function(tracks, minYday, maxYday) {
        tracksDF <- ld(tracks$tracks)
        date <- as.POSIXlt(tracksDF$date)
        years <- date$year + 1900
        randomizedDayTracksDF <- data.frame()
        
        for (year in sort(unique(years))) {
          yearIndex <- years == year
          days <- date$yday[yearIndex]
          daysPeriod <- subset(days, days >= minYday & days <= maxYday)
          randomDay <- sample(daysPeriod, 1)
          message("Randomize day ", randomDay + 1, " for year ", year)
          dayIndex <- days == randomDay
          randomizedDayTracksDF <- rbind(randomizedDayTracksDF, tracksDF[yearIndex,][dayIndex,])
        }

        tracks$tracks <- dl(randomizedDayTracksDF)        
        return(tracks)
      }, minYday=minYday, maxYday=maxYday)
    },
    
    findIntersections.NEW = function(surveyRoutes, dimension, save=FALSE) {
      for (iteration in 1:length(tracksFiles)) {
        tracks <- SimulatedTracks$new(study=study, iteration=iteration)
        tracks$loadTracks()
        
        message("Finding intersections for response = ", study$response, ", iteration = ", tracks$iteration, ", thin = ", tracks$thinId)
        intersections <- SimulatedIntersections$new(study=study, iteration=tracks$iteration)
        intersections$findIntersections(tracks, surveyRoutes, dimension)
        
      }
    },

    findIntersections = function(surveyRoutes, dimension, save=FALSE) {
      intersectionsList <- applyTracks(function(tracks, surveyTracks, dimension, save) {
        message("Finding intersections for response = ", study$response, ", iteration = ", tracks$iteration, ", thin = ", tracks$thinId)
        intersections <- SimulatedIntersections$new(study=study, iteration=tracks$iteration)
        intersections$findIntersections(tracks, surveyRoutes, dimension)
        if (tracks$thinId != 1 & save)
          stop("Saving thinned tracks intersections unsupported.")
        if (save)
          intersections$saveIntersections()
        return(intersections)
      }, surveyTracks=surveyTracks, dimension=dimension, save=save)

      intersectionsCollection <- SimulatedIntersectionsCollection$new(study=study, intersectionsList=intersectionsList)
      return(intersectionsCollection)
    }#,
    
    # Divide-and-conquer approach to avoid taking up all memory
#    thinTracksAndCombineIntersections = function(surveyRoutes, dimension=1, save=FALSE) {
#      sampleIntervals <- MovementSampleIntervals$new(study=study)
#      nIterations <- getNumberOfIterations()
#      intersections <- NULL
#      for (i in 1:nIterations) {
#        message("Processing iteration ", i, " / ", nIterations, "...")
#        thinnedTracks <- sampleIntervals$getThinnedTracksSampleIntervals(getTracks(i))
#        thinnedTracks$getDistances()  
#        x <- thinnedTracks$findIntersections(surveyRoutes=surveyRoutes, dimension=dimension, save=FALSE)  
#        if (is.null(intersections)) intersections <- x else intersections$combineIntersections(x)
#      }
#      if (save) intersections$saveIntersections()
#      return(invisible(intersections))
#    }
  )
)
