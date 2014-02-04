library(sp)

Intersections <- setRefClass(
  Class = "Intersections",
  fields = list(
    study = "Study",
    intersections = "SpatialPointsDataFrame"
  ),
  methods = list(
    initialize = function(preprocessData=FALSE, ...) {
      callSuper(...)
      if (preprocessData) saveIntersections()
      return(.self)
    },
    
    newInstance = function() {
      loadIntersections()
    },
    
    getIntersectionsFileName = function() {
      return(study$context$getFileName(context$resultDataDirectory, name="Intersections", region=study$studyArea$region))
    },
    
    saveIntersections = function() {
      stop("Override saveIntersections() method.")
    },
        
    loadIntersections = function() {
      load(getIntersectionsFileName(), envir=as.environment(.self))
    },
    
    getIntersectionsLocations = function() {
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
    initialize = function(preprocessData=FALSE, ...) {
      if (preprocessData) stop("Unsupported.")
      callSuper(...)
      return(.self)
    },
    
    newInstance = function() {
      callSuper()
    },
    
    findIntersections = function(tracks, surveyRoutes, ...) {
      tracksSP <- tracks$getSpatialLines()      
      findIntersectionsMatrix(tracksSP, surveyRoutes$surveyRoutes, ...)      
      aggregateIntersectionsMatrix(tracks, surveyRoutes)
    },
        
    findIntersectionsMatrix = function(tracks, surveyRoutes, dimension=1) {
      library(sp)
      library(CNPCluster)

      nSurveyRoutes <- length(surveyRoutes)
      nTracks <- length(tracks)
      #intersectionsMatrix <- matrix(0, nrow=nSurveyRoutes, ncol=nTracks)
      
      # Assumes all survey routes were sampled each year
      # TODO: If only one agent, do we still get matrix not vector?
      if (dimension == 1) {
        intersectionsMatrix <<- cnpClusterListApplyGeneric(1:nTracks, function(j, surveyRoutes, tracks, cluster) {          
          library(plyr)
          library(rgeos)
          
          nSurveyRoutes <- length(surveyRoutes)
          nTracks <- length(tracks)
          message("Finding intersections for track ", j," / ", nTracks, "...")
          
          x <- laply(1:nSurveyRoutes,
            function(i, surveyRoutes, tracks, j) {
             return(length(gIntersection(surveyRoutes[i], tracks[j,], byid=TRUE)))
            }, surveyRoutes=surveyRoutes, tracks=tracks, j=j, .parallel=TRUE)
          
          return(as.matrix(x))
        }, surveyRoutes=surveyRoutes, tracks=tracks)
      }
      else if (dimension == 2) {
        intersectionsMatrix <<- cnpClusterListApplyGeneric(1:nSurveyRoutes, function(i, surveyRoutes, tracks) {
          library(plyr)
          library(rgeos)
          
          nSurveyRoutes <- length(surveyRoutes)
          nTracks <- length(tracks)
          message("Finding intersections for survey route ", i, " / ", nSurveyRoutes, "...")
          
          x <- laply(1:nTracks,
            function(j, surveyRoutes, tracks, i) {
             return(length(gIntersection(surveyRoutes[i], tracks[j,], byid=TRUE)))
            }, surveyRoutes=surveyRoutes, tracks=tracks, i=i, .parallel=TRUE)
          
          return(as.matrix(x))
        }, surveyRoutes=surveyRoutes, tracks=tracks)
        
        intersectionsMatrix <<- t(intersectionsMatrix)
      }

      message("Found ", sum(intersectionsMatrix) / nTracks, " intersections per track.")
      message("Found ", sum(intersectionsMatrix) / nSurveyRoutes, " intersections per survey route.")
      
      rownames(intersectionsMatrix) <<- names(surveyRoutes)
      colnames(intersectionsMatrix) <<- sapply(tracks@lines, function(x) x@ID)     
    },
    
    aggregateIntersectionsMatrix = function(tracks, surveyRoutes) {
      tracksDF <- ld(tracks$tracks)
      burstYear <- ddply(tracksDF, .(burst), function(x) {
        date <- as.POSIXlt(x$date)
        data.frame(burst=x$burst[1],
                   year=date$year[1] + 1900,
                   duration=round(as.numeric(difftime(max(x$date), min(x$date))))) # Here it is assumed that the observation period is continuous
      })
      
      data <- data.frame()
      for (year in sort(unique(burstYear$year))) {
        yearToBurstsIndex <- burstYear$year == year
        bursts <- burstYear[yearToBurstsIndex,]$burst
        duration <- burstYear[yearToBurstsIndex,]$duration[1]
        centroids <- coordinates(surveyRoutes$centroids)
        x <- data.frame(x=centroids[,1], y=centroids[,2],
                        year=year,
                        response=study$response,
                        intersections=rowSums(intersectionsMatrix[,bursts,drop=F]),
                        duration=duration,
                        length=surveyRoutes$lengths,
                        distance=tracks$distances)
        data <- rbind(data, x)
      }
      
      intersections <<- SpatialPointsDataFrame(coords=data[,c("x","y")], data=data[,!names(data) %in% c("x","y")], proj4string=surveyRoutes$surveyRoutes@proj4string, match.ID=FALSE)
    },
    
    getIntersectionsFileName = function() {
      if (inherits(study, "undefinedField") | length(iteration) == 0)
        stop("Provide study and iteration parameters.")
      return(study$context$getLongFileName(dir=study$context$resultDataDirectory, name="Intersections", response=study$response, region=study$studyArea$region, tag=iteration))
    },
    
    saveIntersections = function() {
      fileName <- getIntersectionsFileName()
      message("Saving intersections to ", fileName)
      save(intersections, intersectionsMatrix, iteration, file=fileName)
    }
  )
)

FinlandWTCIntersections <- setRefClass(
  Class = "FinlandWTCIntersections",
  contains = "Intersections",
  methods = list(
    newInstance = function(...) {
      callSuper(...)
      #studyArea <<- FinlandStudyArea(context=context)$newInstance()
      return(.self)
    },
    
    saveIntersections = function() {
      library(gdata)
      library(sp)
      
      wtc <- read.xls(file.path(study$context$rawDataDirectory, "Kolmiot_1989_2011.xls"))
      names(wtc)[1:9] <- c("id","y","x","length","forest","field","year","date","duration")
      names(wtc)[16] <- "canis.lupus"
      names(wtc)[29] <- "lynx.lynx"
      names(wtc)[33] <- "rangifer.tarandus.fennicus"
      wtc$date <- strptime(wtc$date, "%Y%m%d")
      wtc$x <- (wtc$x+3000)*1000
      wtc$y <- wtc$y*1000
      wtc <- subset(wtc, duration<4 & duration>0 & length>0)  
      # Transects 168 and 1496 are at the same location but separate, take average, ignore or something...
      wtc <- subset(wtc, id!=1496)
      wtc <- SpatialPointsDataFrame(cbind(wtc$x, wtc$y), data=wtc, proj4string=CRS("+init=epsg:2393"))
      
      intersections <<- wtc[,1:9]
      intersections <<- rbind(intersections, intersections, intersections)
      intersections$response <<- c(rep("canis.lupus", times=nrow(wtc)), rep("lynx.lynx", times=nrow(wtc)), rep("rangifer.tarandus.fennicus", times=nrow(wtc)))
      intersections$intersections <<- c(wtc$canis.lupus, wtc$lynx.lynx, wtc$rangifer.tarandus.fennicus)
      
      save(intersections, file=getDataFileName())
    }
  )
)
