library(adehabitatLT)
setOldClass("ltraj")

Tracks <- setRefClass(
  Class = "Tracks",
  fields = list(
    study = "Study",
    tracks = "ltraj",
    distance = "numeric"
  ),
  methods = list(
    initialize = function(preprocessData=FALSE, ...) {
      callSuper(...)
      if (preprocessData) saveTracks()
      return(.self)
    },
    
    getTracksDirectory = function() return(study$context$scratchDirectory),
    
    getTracksFileName = function() {
      return(study$context$getFileName(getTracksDirectory(), name="Tracks", response=study$response, region=study$studyArea$region))
    },
    
    saveTracks = function() {
      stop("Override saveData() method.")
    },
    
    loadTracks = function() {
      library(adehabitatLT)
      load(getTracksFileName(), envir=as.environment(.self))
      return(invisible(.self))
    },
    
    getSpatialLines = function() {
      library(sp)
      tracksSP <- ltraj2sldf(tracks, byid=FALSE)
      proj4string(tracksSP) <- study$studyArea$proj4string
      return(tracksSP)
    },
    
    plotTracks = function(surveyRoutes) {
      library(sp)
      library(adehabitatLT)
      
      tracksDF <- ld(tracks)
      plot(tracksDF$x, tracksDF$y, type="n")
      apply(fun=function(x) lines(x$x, x$y, col=x$id), variables=.(burst))
      plot(study$studyArea$boundary, add=T)
      if (!missing(surveyRoutes))
        plot(surveyRoutes$surveyRoutes, col="blue", add=T)
    },
    
    apply = function(variables=.(id), fun, ..., combine=F) {
      tracksDF <- ld(tracks)
      result <- dlply(.data=tracksDF, .variables=variables, .fun=fun, ...)
      if (combine) result <- do.call("rbind", result)
      return(result)
    },
        
    getSampleIntervals = function() {
      intervals <- MovementSampleIntervals$new()
      intervals$getSampleIntervals(tracks=.self)
      return(intervals)
    },
    
    getThinnedTracksSampleIntervals = function() {
      intervals <- TracksSampleIntervals$new()
      intervals$getThinnedTracksSampleIntervals(tracks=.self)
      return(intervals)
    },
    
    # Note! Call this before randomizing observation days. Otherwise you'll lose details of the last movements.
    determineDistances = function() {
      warning("Spatial variation in distances not considered in vicinity of the survey routes.")
      
      tracksDF <- ld(tracks)
      d <- as.POSIXlt(tracksDF$date)
      tracksDF$yday <- d$yday
      tracksDF$year <- d$year
      distances <- ddply(.data=tracksDF, .variables=.(id, burst, yday, year), .fun=function(x) {
        s <- sum(x$dt, na.rm=T) / 3600
        if (s < 23 | s > 25) return(NA)
        return(sum(x$dist) / sum(x$dt) * 24 * 3600)
      })$V1
      
      if (all(is.na(distances)))
        stop("Unable to determine movement distance.")
      
      message(sum(!is.na(distances)), " / ", length(distances), " of day movements used to determine mean movement distance = ", mean(distances, na.rm=T), " ± ", sd(distances, na.rm=T), ".")
      distance <<- mean(distances, na.rm=T)
      #distances <<- rep(distance, times=length(surveyRoutes$surveyRoutes))
      return(invisible(distances))
    },
    
    getHabitatPreferences = function(habitatWeightsTemplate, nSamples=30, save=FALSE) {
      habitatPreferences <- HabitatSelection$new(study=study)
      habitatPreferences$getHabitatPreferences(tracks=.self, habitatWeightsTemplate=habitatWeightsTemplate, nSamples=nSamples, save=save)
      return(habitatPreferences)
    }
  )
)

SimulatedTracks <- setRefClass(
  Class = "SimulatedTracks",
  contains = "Tracks",
  fields = list(
    iteration = "integer",
    thinId = "integer"
  ),
  methods = list(
    initialize = function(xy, id, date, thinId, preprocessData=FALSE, ...) {
      if (!missing(xy) && !missing(id) && !missing(date)) setTracks(xy=xy, id=id, date=date)
      if (missing(thinId)) thinId <<- as.integer(1)
      else thinId <<- thinId
      callSuper(preprocessData=preprocessData, ...)
    },
    
    setTracks = function(xy, id, date) {
      library(adehabitatLT)
      d <- as.POSIXlt(date)
      burst <- paste(id, d$year, sep=".")
      tracks <<- as.ltraj(xy=xy, date=date, id=id, burst=burst)
    },
    
    getTracksFileName = function() {
      if (inherits(study, "undefinedField") | length(iteration) == 0)
        stop("Provide study and iteration parameters.")
      return(study$context$getLongFileName(dir=getTracksDirectory(), name="Tracks", response=study$response, region=study$studyArea$region, tag=iteration))
    },
    
    saveTracks = function() {
      fileName <- getTracksFileName()
      message("Saving tracks to ", fileName)
      if (thinId != 1)
        stop("Saving thinned tracks unsupported.")
      save(tracks, iteration, thinId, file=fileName)
    },

    .internal.thin = function(by) {
      library(plyr)
      
      tracksDF <- ld(tracks)
      d <- as.POSIXlt(tracksDF$date)
      tracksDF$yday <- d$yday
      tracksDF$year <- d$year
      
      oldDt <- mean(tracksDF$dt, na.rm=T)
      
      tracksDFThinned <- ddply(tracksDF, .variables=.(id, burst, year, yday), .fun=function(x, by) {
        n <- nrow(x)
        if (n==1) return(NULL)
        retainIndex <- seq(1, n, by=by) 
        return(x[retainIndex,])
      }, by=by)
      
      if (nrow(tracksDFThinned) == 0) return(NULL)
      
      tracksThinned <- dl(tracksDFThinned)
      tracksDFThinned <- ld(tracksThinned)
      
      newDt <- mean(tracksDFThinned$dt, na.rm=T)
      message("Thinned movements from dt = ", oldDt, " to dt = ", newDt, " (note: these values can be misleading)")      
      return(tracksThinned)
    },
    
    thin = function(by) {
      message("Thinning iteration = ", iteration, ", thin = ", thinId)
      thinnedTracksDF <- .internal.thin(by=by)
      if (is.null(thinnedTracksDF)) return(NULL)
      thinnedTracks <- SimulatedTracks$new(study=study, tracks=thinnedTracksDF, iteration=iteration, thinId=as.integer(thinId+1))
      return(thinnedTracks)
    }
  )
)

FinlandWTCTracks <- setRefClass(
  Class = "FinlandWTCTracks",
  contains = "Tracks",
  fields = list(
  ),
  methods = list(
    saveTracks = function() {
      maxSpeed <- 40 # km/h
      
      library(sp)
      library(plyr)
      library(rgdal)
      library(adehabitatLT)
      
      gps <- if (study$response == "canis.lupus") preprocessWolfData()
      else if (study$response == "lynx.lynx") preprocessLynxData()
      else if (study$response == "rangifer.tarandus.fennicus") preprocessReindeerData()
      
      fromProj <- gps$fromProj
      gps <- gps$gps
      
      message("Cleaning data...")
      gps <- ddply(gps, .(id), function(x) x[!duplicated(x$date),]) # Remove duplicate dates
      gps <- gps[complete.cases(gps),] # Remove missing data
      coordinates(gps) <- ~x+y
      proj4string(gps) <- fromProj
      gps <- spTransform(gps, study$studyArea$proj4string)
      gpsLT <- as.ltraj(coordinates(gps), gps$date, gps$id, infolocs=NULL)
      gpsLT <- gdltraj(gpsLT, min=0, max=2, type="mon") # Limit to WTC period, TODO: period is 3 months in the Northern Finland
      gpsLT <- cutltraj(gpsLT, "dt > 3600*24*200", nextr=TRUE) # Cut years apart
      gpsLT <- cutltraj(gpsLT, paste("(dist/1000)/(dt/3600) >", maxSpeed)) # Remove unrealistic movements
      assign(".b", study$studyArea$boundary@polygons[[1]]@Polygons[[1]]@coords, envir=.GlobalEnv)
      gpsLT <- cutltraj(gpsLT, "point.in.polygon(x, y, .b[,1], .b[,2])==0") # Remove movements outside the boundary
      rm(".b", envir=.GlobalEnv)
      
      message("Saving processed data...")
      tracks <<- gpsLT
      save(tracks, file=getTracksFileName())
      return(invisible(.self))
    },
    
    preprocessWolfData = function() {
      library(gdata)
      
      message("Reading raw data...")
      gps.raw <- read.xls(file.path(study$context$rawDataDirectory, "GPS_Finland_GPS_Finland_RKTL_Wolf_ec_import_tracks.xlsx"), na.strings="")  
      gps <- data.frame(id=substr(gps.raw$UnitID, 11, length(gps.raw$UnitID)),
                        date=as.POSIXct(gps.raw$Time*24*60*60, origin=as.POSIXct(gps.raw$Date, format="%Y-%m-%d")),
                        x=as.numeric(as.character(gps.raw$X)), y=as.numeric(as.character(gps.raw$Y)))
      index <- gps$x < 0 | gps$x > 200 | gps$y < 0 | gps$y > 200
      gps <- gps[!index,]
      
      return(list(gps=gps, fromProj="+init=epsg:4326"))
    },
    
    preprocessLynxData = function() {
      library(gdata)
      
      message("Reading raw data...")
      gps.raw <- read.xls(file.path(study$context$rawDataDirectory, "GPS_Finland_GPS_Finland_RKTL_Lynx_ec_import_tracks.xls"))  
      gps <- data.frame(id=gps.raw$UnitID,
                        date=as.POSIXct(paste(gps.raw$Date, gps.raw$Time), format="%m/%d/%Y %I:%M:%S %p"),
                        x=gps.raw$X, y=gps.raw$Y)
      
      return(list(gps=gps, fromProj="+init=epsg:2393"))
    },
    
    preprocessReindeerData = function() {
      library(gdata)
      
      message("Reading raw data...")
      gps.raw <- read.xls(file.path(study$context$rawDataDirectory, "GPS_Finland_GPS_Finland_RKTL_ForestReindeer_ec_import_tracks.xlsx"))
      gps <- data.frame(id=gps.raw$UnitID,
                        date=as.POSIXct(gps.raw$Date*60*60*24-60*60*2*24-60*60*.3, origin="1900-01-01"),
                        x=gps.raw$X, y=gps.raw$Y)
      
      return(list(gps=gps, fromProj="+init=epsg:4326"))
    }
  )
)
