library(sp)

Intersections <- setRefClass(
  Class = "Intersections",
  fields = list(
    context = "Context",
    studyArea = "StudyArea",
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
      return(context$getFileName(context$processedDataDirectory, name="WTC", region=studyArea$region))
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
  methods = list(
    initialize = function(preprocessData=FALSE, ...) {
      callSuper(...)
      if (preprocessData) saveIntersections()
      return(.self)
    },
    
    newInstance = function() {
      callSuper()
    },
    
    findIntersections = function(tracks, surveyRoutes, ...) {
      tracksSP <- tracks$getSpatialLines()
      intersectionsMatrix <- findIntersections.Internal(tracksSP, surveyRoutes$surveyRoutes, ...)
      
      #tracksDF <- ld(tracks$tracks)
      #ddply(tracksDF, .(burst), function(x) {
      #  date <- as.POSIXlt(x$date)
      #  data.frame(year=date$year)
      #}
      #intersections <<- data.frame()
      
      return(intersectionsMatrix)
    },
    
    findIntersections.Internal = function(tracks, surveyRoutes, dimension=1) {
      library(sp)
      library(CNPCluster)
      
      nSurveyRoutes <- length(surveyRoutes)
      nTracks <- length(tracks)
      #intersectionsMatrix <- matrix(0, nrow=nSurveyRoutes, ncol=nTracks)
      
      if (dimension == 1) {
        intersectionsMatrix <- cnpClusterListApplyGeneric(1:nTracks, function(j, surveyRoutes, tracks, cluster) {          
          library(plyr)
          library(rgeos)
          
          nSurveyRoutes <- length(surveyRoutes)
          nTracks <- length(tracks)
          message("Finding intersections for track ", j," / ", nTracks, "...")
          
          x <- laply(1:nSurveyRoutes,
            function(i, surveyRoutes, tracks, j) {
             return(length(gIntersection(surveyRoutes[i], tracks[j,], byid=TRUE)))
            }, surveyRoutes=surveyRoutes, tracks=tracks, j=j, .parallel=TRUE)
          
          return(x)
        }, surveyRoutes=surveyRoutes, tracks=tracks)
      }
      else if (dimension == 2) {
        intersectionsMatrix <- cnpClusterListApplyGeneric(1:nSurveyRoutes, function(i, surveyRoutes, tracks) {
          library(plyr)
          library(rgeos)
          
          nSurveyRoutes <- length(surveyRoutes)
          nTracks <- length(tracks)
          message("Finding intersections for survey route ", i, " / ", nSurveyRoutes, "...")
          
          x <- laply(1:nTracks,
            function(j, surveyRoutes, tracks, i) {
             return(length(gIntersection(surveyRoutes[i], tracks[j,], byid=TRUE)))
            }, surveyRoutes=surveyRoutes, tracks=tracks, i=i, .parallel=TRUE)
          
          return(x)
        }, surveyRoutes=surveyRoutes, tracks=tracks)
        
        intersectionsMatrix <- t(intersectionsMatrix)
      }

      rownames(intersectionsMatrix) <- names(surveyRoutes)
      colnames(intersectionsMatrix) <- sapply(tracks@lines, function(x) x@ID)
     
      message("Found ", sum(intersectionsMatrix) / nTracks, " intersections per track.")
      message("Found ", sum(intersectionsMatrix) / nSurveyRoutes, " intersections per survey route.")
            
      return(invisible(intersectionsMatrix))
    }
  )
)

FinlandWinterTrackCounts <- setRefClass(
  Class = "FinlandWinterTrackCounts",
  contains = "Intersections",
  methods = list(
    newInstance = function() {
      callSuper()
      studyArea <<- FinlandStudyArea(context=context)$newInstance()
      return(.self)
    },
    
    saveIntersections = function() {
      library(gdata)
      library(sp)
      
      wtc <- read.xls(file.path(context$rawDataDirectory, "Kolmiot_1989_2011.xls"))
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
      intersections <<- wtc
      save(intersections, file=getDataFileName())
    }
  )
)
