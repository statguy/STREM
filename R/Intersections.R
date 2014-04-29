library(sp)

Intersections <- setRefClass(
  Class = "Intersections",
  fields = list(
    study = "Study",
    intersections = "SpatialPointsDataFrame"
  ),
  methods = list(
    initialize = function(...) {
      callSuper(...)
      return(invisible(.self))
    },
    
    getCoordinates = function() return(coordinates(intersections)),
    getData = function() return(intersections@data),
    
    getIntersectionsFileName = function(tag) {
      if (inherits(study, "undefinedField"))
        stop("Provide study parameters.")
      if (missing(tag))
        return(study$context$getFileName(dir=study$context$resultDataDirectory, name="Intersections", response=study$response, region=study$studyArea$region))  
      else
        return(study$context$getLongFileName(dir=study$context$resultDataDirectory, name="Intersections", response=study$response, region=study$studyArea$region, tag=tag))
    },
    
    saveIntersections = function(fileName) {
      stop("Override saveIntersections() method.")
    },
    
    loadIntersections = function(fileName=getIntersectionsFileName()) {
      load(fileName, envir=as.environment(.self))
      return(invisible(.self))
    },
    
    getSurveyLocations = function() {
      return(SpatialPoints(unique(coordinates(intersections)), proj4string=intersections@proj4string))
    }
  )
)

SimulatedIntersections <- setRefClass(
  Class = "SimulatedIntersections",
  contains = "Intersections",
  fields = list(
    intersectionsMatrix = "matrix",
    iteration = "integer"
  ),
  methods = list(
    initialize = function(...) {
      callSuper(...)
      return(invisible(.self))
    },
    
    findIntersections = function(tracks, surveyRoutes, ...) {
      tracksSP <- tracks$getSpatialLines()
      message("Counting intersections...")
      findIntersectionsMatrix(tracksSP, surveyRoutes$surveyRoutes, ...)
      nSurveyRoutes <- length(surveyRoutes$surveyRoutes)
      nTracks <- length(tracksSP)
      nYears <- length(unique(tracks$tracks$year))
      message("Found ", sum(intersectionsMatrix) / nTracks / nYears, " intersections per track per year.")
      message("Found ", sum(intersectionsMatrix) / nSurveyRoutes / nYears, " intersections per survey route per year.")
      aggregateIntersectionsMatrix(tracks, surveyRoutes)
    },
    
    # TODO: Support for WTC survey routes for each year
    findIntersectionsMatrix = function(tracks, surveyRoutes, dimension=1) {
      library(sp)
      
      nSurveyRoutes <- length(surveyRoutes)
      nTracks <- length(tracks)
      #intersectionsMatrix <<- matrix(0, nrow=nSurveyRoutes, ncol=nTracks)
      
      if (dimension == 1) {
        # Assumes that the survey routes are the same every year
        countIntersections <- function(i, surveyRoutes, tracks) {
          library(plyr)
          #library(rgeos)

          nSurveyRoutes <- length(surveyRoutes)
          nTracks <- length(tracks)
          message("Finding intersections for survey route ", i, " / ", nSurveyRoutes, " for ", nTracks, " tracks...")
                   
#           i <- 2
#           count <- 0
#           for (j in 1:length(tracks)) {
#             if (nrow(tracks[j]@lines[[1]]@Lines[[1]]@coords) == 1) next
#             x <- length(gIntersection(surveyRoutes[i], tracks[j], byid=TRUE))
#             if (x>0) {
#               y <- gIntersection(surveyRoutes[i], tracks[j], byid=TRUE)
#               message(i, ", ", j)
#             }
#             count <- count + x
#           }
#           count
          
          surveyRouteCoords <- coordinates(surveyRoutes[i])[[1]][[1]]
          
          x <- integer(nTracks)
          for (j in 1:nTracks) {
            trackCoords <- coordinates(tracks[j])[[1]][[1]]
            x[j] <- internal_countIntersections(trackCoords, surveyRouteCoords)
          }
          
          return(x)
        }
        
        intersectionsMatrix <<- laply(1:nSurveyRoutes, countIntersections, surveyRoutes=surveyRoutes, tracks=tracks, .drop=FALSE, .parallel=TRUE, .inform=FALSE)
      }
      
      if (F) {
      # Assumes that the survey routes are the same every year
      if (dimension == 1) {
        countIntersections <- function(i, surveyRoutes, tracks) {
          library(plyr)
          library(rgeos)
          
          nSurveyRoutes <- length(surveyRoutes)
          nTracks <- length(tracks)
          message("Finding intersections for survey route ", i, " / ", nSurveyRoutes, " for ", nTracks, " tracks...")
                    
          x <- laply(1:nTracks,
                     function(j, surveyRoutes, tracks, i) {
                       return(length(gIntersection(surveyRoutes[i], tracks[j], byid=TRUE)))
                     }, surveyRoutes=surveyRoutes, tracks=tracks, i=i, .inform=FALSE)
          
          return(x)
        }
        
        intersectionsMatrix <<- laply(1:nSurveyRoutes, countIntersections, surveyRoutes=surveyRoutes, tracks=tracks, .drop=FALSE, .parallel=TRUE, .inform=FALSE)
      }
      else if (dimension == 2) {
        countIntersections <- function(j, surveyRoutes, tracks) {
          library(plyr)
          library(rgeos)
          
          nSurveyRoutes <- length(surveyRoutes)
          nTracks <- length(tracks)
          message("Finding intersections for track ", j," / ", nTracks, " for ", nSurveyRoutes, " survey routes...")
          
          x <- laply(1:nSurveyRoutes,
                     function(i, surveyRoutes, tracks, j) {
                       return(length(gIntersection(surveyRoutes[i], tracks[j], byid=TRUE)))
                     }, surveyRoutes=surveyRoutes, tracks=tracks, j=j, .inform=FALSE)
          
          return(x)
        }
        
        intersectionsMatrix <<- laply(1:nTracks, countIntersections, surveyRoutes=surveyRoutes, tracks=tracks, .drop=FALSE, .parallel=TRUE, .inform=FALSE)
        intersectionsMatrix <<- t(intersectionsMatrix)
      }
      }
      
      rownames(intersectionsMatrix) <<- names(surveyRoutes)
      colnames(intersectionsMatrix) <<- sapply(tracks@lines, function(x) x@ID)     
    },
    
    aggregateIntersectionsMatrix = function(tracks, surveyRoutes) {
      #if (length(tracks$distance) == 0)
      #  stop("Did you forgot to run getDistances() for tracks first?")
      
      tracksObj <- tracks$getTracks()
      tracksDF <- if (inherits(tracksObj, "ltraj")) ld(tracksObj) else tracksObj
      
      burstYear <- ddply(tracksDF, .(id,year), function(x) {
        date <- as.POSIXlt(x$date)
        data.frame(burst=x$burst[1],
                   year=date$year[1] + 1900,
                   duration=round(as.numeric(difftime(max(x$date), min(x$date), units="days")))) # Here it is assumed that the observation period is continuous
      })
      
      data <- data.frame()
      for (year in sort(unique(burstYear$year))) {
        yearToBurstsIndex <- burstYear$year == year
        bursts <- as.character(burstYear[yearToBurstsIndex,]$burst)
        duration <- burstYear[yearToBurstsIndex,]$duration[1]
        centroids <- coordinates(surveyRoutes$centroids)
        x <- data.frame(surveyRoute=rownames(intersectionsMatrix[,bursts,drop=F]),
                        x=centroids[,1], y=centroids[,2],
                        year=year,
                        response=study$response,
                        intersections=rowSums(intersectionsMatrix[,bursts,drop=F]),
                        duration=duration,
                        length=surveyRoutes$lengths,
                        distance=1) # Correct for distance after estimation if not set elsewhere
        data <- rbind(data, x)
      }
      
      intersections <<- SpatialPointsDataFrame(coords=data[,c("x","y")],
                                               data=data[,!names(data) %in% c("x","y")],
                                               proj4string=surveyRoutes$surveyRoutes@proj4string,
                                               match.ID=FALSE)
    },
    
    getIntersectionsFileName = function(tag) {
      if (inherits(study, "undefinedField") | length(iteration) == 0)
        stop("Provide study and iteration parameters.")
      if (missing(tag))
        return(study$context$getLongFileName(dir=study$context$resultDataDirectory, name="Intersections", response=study$response, region=study$studyArea$region, tag=iteration))
      else
        return(study$context$getLongFileName(dir=study$context$resultDataDirectory, name="Intersections", response=study$response, region=study$studyArea$region, tag=paste(iteration, tag, sep="-")))
    },
    
    saveIntersections = function(fileName=getIntersectionsFileName()) {
      message("Saving intersections to ", fileName)
      save(intersections, intersectionsMatrix, iteration, file=fileName)
    }
  )
)

FinlandWTCIntersections <- setRefClass(
  Class = "FinlandWTCIntersections",
  contains = c("Intersections", "FinlandCovariates"),
  fields = list(
    maxDuration = "numeric"
  ),
  methods = list(
    initialize = function(maxDuration=Inf, ...) {
      callSuper(covariatesName="FinlandWTCIntersectionsCovariates", ...)
      maxDuration <<- maxDuration
      return(invisible(.self))
    },
    
    saveIntersections = function() {
      library(gdata)
      library(sp)
      
      message("Processing ", study$response, "...")
      
      wtc <- read.xls(file.path(study$context$rawDataDirectory, "Kolmiot_1989_2011.xls"))
      names(wtc)[1:9] <- c("id","y","x","length","forest","field","year","date","duration")
      names(wtc)[16] <- "canis.lupus"
      names(wtc)[29] <- "lynx.lynx"
      names(wtc)[33] <- "rangifer.tarandus.fennicus"
      
      responseIndex <- switch (study$response,
              canis.lupus=16,
              lynx.lynx=29,
              rangifer.tarandus.fennicus=33)
      wtc <- wtc[,c(1:9,responseIndex)]
      colnames(wtc)[ncol(wtc)] <- "intersections"
      wtc$response <- study$response
      wtc <- subset(wtc, response==response)
      
      wtc$date <- strptime(wtc$date, "%Y%m%d")
      date <- as.POSIXlt(wtc$date)
      index <- date$year == 0
      date$year[index] <- 100
      wtc$date <- as.POSIXct(date)
      
      wtc$x <- (wtc$x+3000)*1000
      wtc$y <- wtc$y*1000
      wtc <- subset(wtc, duration<4 & duration>0 & length>0)  
      # Transects 168 and 1496 are at the same location but separate, take average, ignore or something...
      wtc <- subset(wtc, id!=1496)
      xy <- cbind(wtc$x, wtc$y)
      wtc$x <- NULL
      wtc$y <- NULL
      wtc$length <- wtc$length * 1000
      
      b <- nrow(wtc)
      nonZero <- wtc[wtc$duration > maxDuration,]$intersections
      nonZero <- sum(nonZero != 0) / length(nonZero)
      retainIndex <- wtc$duration <= maxDuration
      wtc <- wtc[retainIndex,]
      a <- nrow(wtc)
      message("Removed ", round((1-a/b)*100), "% of intersections with duration > ",
              maxDuration, " of which ", round(nonZero*100), "% are non-zero.")
      
      intersections <<- SpatialPointsDataFrame(xy[retainIndex,], data=wtc, proj4string=CRS("+init=epsg:2393"))

      save(intersections, file=getIntersectionsFileName())
    },
    
    predictDistances = function(formula=study$getDistanceCovariatesModel(), intervalH=study$getTrackSampleInterval()) {
      intersections$distance <<- study$predictDistances(formula=formula, data=covariates, intervalH=intervalH)      
      return(invisible(.self))
    }
  )
)
