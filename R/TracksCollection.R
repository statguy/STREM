TracksCollection <- setRefClass(
  Class = "TracksCollection",
  fields = list(
    study = "Study",
    tracksList = "list"
  ),
  methods = list(
    initialize = function(...) {
      callSuper(...)
    },
  
    length = function() return(base::length(tracksList)),
    
    add = function(tracks) {
      tracksList[[length() + 1]] <<- tracks
    },
    
    get = function(index) {
      if (base::length(tracksList) < index) stop("Invalid index = ", index, ".")
      return(tracksList[[index]])
    },
    
    load = function() {
      stop("Unimplemented method.")
    },
    
    apply = function(fun, ...) {
      return(lapply(X=tracksList, FUN=fun, ...))
    },
    
    findIntersections = function(surveyRoutes, dimension) {
      stop("Unimplemented method.")
    }
  )
)

SimulatedTracksCollection <- setRefClass(
  Class = "SimulatedTracksCollection",
  contains = "TracksCollection",
  methods = list(
    load = function() {     
      tracksFiles <- study$context$listFiles(dir=study$context$processedDataDirectory, name="Tracks", response=paste(study$response, "\\d+", sep="-"), region=study$studyArea$region)
      
      for (iteration in 1:base::length(tracksFiles)) { # TODO: better to use iteration number rather than number of files
        tracks <- SimulatedTracks$new(study=study, iteration=iteration)
        tracks$loadTracks()
        add(tracks)
      }
    },
    
    randomizeObservationDays = function(minYday=0, maxYday=60) {
      tracks <- apply(function(tracks, minYday, maxYday) {
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
    
    findIntersections = function(surveyRoutes, dimension, save=FALSE) {
      intersectionsList <- apply(function(tracks, surveyTracks, dimension, save) {
        message("Finding intersections for response = ", study$response, ", iteration = ", tracks$iteration)        
        intersections <- SimulatedIntersections$new(study=study, iteration=tracks$iteration)
        intersections$findIntersections(tracks, surveyRoutes, dimension)
        if (save) intersections$saveIntersections()
        return(intersections)
      }, surveyTracks=surveyTracks, dimension=dimension, save=save)
      
      intersectionsCollection <- IntersectionsCollection$new(study=study, intersectionsList=intersectionsList)
      return(intersectionsCollection)
    }
  )
)
