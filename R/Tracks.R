library(adehabitatLT)
setOldClass("ltraj")

Tracks <- setRefClass(
  Class = "Tracks",
  fields = list(
    context = "Context",
    study = "Study",
    tracks = "ltraj"
  ),
  methods = list(
    initialize = function(preprocessData=FALSE, ...) {
      callSuper(...)
      if (preprocessData) saveTracks()
      return(.self)
    },
    
    getDataFileName = function() {
      return(context$getFileName(context$processedDataDirectory, name="Tracks", response=study$response, region=study$studyArea$region))
    },
    
    saveTracks = function() {
      stop("Override saveData() method.")
    },
    
    loadTracks = function() {
      load(getDataFileName(), envir=as.environment(.self))
      return(tracks)
    },
    
    plotTracks = function() {
      library(sp)
      library(adehabitatLT)
      
      tracksDF <- ld(tracks)
      plot(tracksDF$x, tracksDF$y, type="n")
      apply(fun=function(x) lines(x$x, x$y, col=x$id))
      plot(study$studyArea$boundary, add=T)
    },
    
    apply = function(variables=.(id), fun, ..., combine=F) {
      tracksDF <- ld(tracks)
      result <- dlply(.data=tracksDF, .variables=variables, .fun=fun, ...)
      if (combine) result <- do.call("rbind", result)
      return(result)
    },
        
    thin = function(by) {
      library(plyr)
      
      tracksDF <- ld(tracks)
      d <- as.POSIXlt(tracksDF$date)
      tracksDF$yday <- d$yday
      tracksDF$year <- d$year
      
      tracksDFThinned <- ddply(tracksDF, .variables=.(id, year, yday), .fun=function(x, by) {
        n <- nrow(x)
        if (n==1) return(NULL)
        i <- seq(1, n, by=by)
        return(x[i,])
      }, by=by)
      
      if (nrow(tracksDFThinned) == 0) return(NULL)
      tracksThinned <- dl(tracksDFThinned)
      return(Tracks$new(context=context, study=study, tracks=tracksThinned))
    },
    
    getSampleIntervals = function() {
      intervals <- TracksSampleIntervals$new()
      intervals$getSampleIntervals(tracks=.self)
      return(intervals)
    },
    
    getThinnedTrackSampleIntervals = function(maxThins) {
      intervals <- TracksSampleIntervals$new()
      intervals$getThinnedTracksSampleIntervals(tracks=.self, maxThins=maxThins)
      return(intervals)
    }    
  )
)

SimulatedTracksCollection <- setRefClass(
  Class = "TracksCollection",
  fields = list(
    tracksList = "list"
  ),
  methods = list(
    initialize = function() {
      callSuper(tracksList=list())
    },
    
    length = function() return(base::length(tracksList)),
    
    add = function(tracks) {
      tracksList[[length() + 1]] <<- tracks
    },
    
    get = function(index) {
      if (base::length(tracksList) < index) stop("Invalid index = ", index, ".")
      return(tracksList[[index]])
    },
    
    load = function(study) {     
      context <- study$context
      tracksFiles <- context$listFiles(dir=context$processedDataDirectory, name="Tracks", response=paste(study$response, "\\d+", sep="-"), region=study$studyArea$region)
      
      for (iteration in 1:base::length(tracksFiles)) { # TODO: better to use iteration number rather than number of files
        tracks <- SimulatedTracks$new(context=study$context, study=study, iteration=iteration)
        tracks$loadTracks()
        add(tracks)
      }
    },
    
    apply = function(fun, ...) {
      return(lapply(X=tracksList, FUN=fun, ...))
    }
  )
)

TracksSampleIntervals <- setRefClass(
  Class = "TracksSampleIntervals",
  fields = list(
    intervals = "data.frame"
  ),
  methods = list(
    initialize = function() {
      callSuper(intervals=data.frame())
    },
    
    # Determines sampling intervals of the recorded movements for each day
    getSampleIntervals = function(tracks) {
      library(plyr)
      
      tracksDF <- ld(tracks$tracks)
      tracksDF$row <- 1:nrow(tracksDF)
      date <- as.POSIXlt(tracksDF$date)
      tracksDF$yday <- date$yday
      tracksDF$year <- date$year
      intervals <<- ddply(tracksDF, .(burst, yday, year), function(x) {
        s <- sum(x$dt, na.rm=T) / 3600
        if (s < 23 | s > 25) return(NULL)
        distKm <- sum(x$dist) / 1e3
        if (is.na(distKm) | distKm > 100) return(NULL)
        intervalMin <- 24 / nrow(x) * 60
        if (intervalMin > 24*60) return(NULL)
        x$intervalMin <- intervalMin
        x$intervalH <- intervalMin / 60
        x$distanceKm <- distKm
        return(x[1, c("id","date","intervalH","intervalMin","distanceKm")])
      })
      
      if (nrow(intervals) == 0) warning("Unable to determine sampling intervals.")
    },
    
    getThinnedTracksSampleIntervals = function(tracks, maxThins) {
      thinnedTracksCollection <- SimulatedTracksCollection$new()
      thinnedTracksCollection$add(tracks)
      intervalsList <- list()
      intervalsList[[1]] <- tracks$getSampleIntervals()
      
      for (i in 2:maxThins) {
        thinnedTracks <- thinnedTracksCollection$get(1)$thin(by=i)
        if (is.null(thinnedTracks)) break
        intervalsList[[i]] <- thinnedTracks$getSampleIntervals()
        thinnedTracksCollection$add(tracks)
      }
      
      intervals <<- ldply(intervalsList, function(x) x$intervals)
      
      return(invisible(thinnedTracksCollection))
    },
    
    plotIntervalDistance = function() {
      library(ggplot2)
      p <- ggplot(intervals, aes(intervalH, distanceKm)) +
        geom_point(aes(color=id)) +
        #geom_line(aes(intervalH, fitted), color="red") +
        ylab("Distance / day (km)") + xlab("Sampling interval (h)") + theme_bw(18)
      print(p)
    }
  )
)

SimulatedTracks <- setRefClass(
  Class = "SimulatedTracks",
  contains = "Tracks",
  fields = list(
    iteration = "integer"
  ),
  methods = list(
    initialize = function(tracksDF, preprocessData=FALSE, ...) {
      if (!missing(tracksDF)) setTracks(tracksDF)
      callSuper(preprocessData=preprocessData, ...)
    },

    setTracks = function(tracksDF) {
      if (inherits(tracksDF, "uninitializedField"))
        stop("Provide tracksDF argument.")
      library(adehabitatLT)      
      xy <- tracksDF[,c("x","y")]
      date <- as.POSIXct(strptime(paste(2000+tracksDF$year, tracksDF$day, tracksDF$hour, "0", "0"), format="%Y %j %H %M %S"))
      id <- tracksDF$agent
      tracks <<- as.ltraj(xy=xy, date=date, id=id, infolocs=NULL)
      iteration <<- tracksDF$iteration[1]      
    },
    
    saveTracks = function() {
      fileName <- getDataFileName()
      message("Saving tracks to ", fileName)
      save(tracks, file=fileName)
    },
    
    getDataFileName = function() {
      if (inherits(study, "undefinedField") | length(iteration) == 0)
        stop("Provide study and iteration parameters.")
      response <- paste(study$response, iteration, sep="-")
      return(context$getFileName(dir=context$processedDataDirectory, name="Tracks", response=response, region=study$studyArea$region))
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
      gps <- spTransform(gps, studyArea$proj4string)
      gpsLT <- as.ltraj(coordinates(gps), gps$date, gps$id, infolocs=NULL)
      gpsLT <- gdltraj(gpsLT, min=0, max=2, type="mon") # Limit to WTC period, TODO: period is 3 months in the Northern Finland
      gpsLT <- cutltraj(gpsLT, "dt > 3600*24*200", nextr=TRUE) # Cut years apart
      gpsLT <- cutltraj(gpsLT, paste("(dist/1000)/(dt/3600) >", maxSpeed)) # Remove unrealistic movements
      assign(".b", studyArea$boundary@polygons[[1]]@Polygons[[1]]@coords, envir=.GlobalEnv)
      gpsLT <- cutltraj(gpsLT, "point.in.polygon(x, y, .b[,1], .b[,2])==0") # Remove movements outside the boundary
      rm(".b", envir=.GlobalEnv)
      
      message("Saving processed data...")
      tracks <<- gpsLT
      save(tracks, file=getDataFileName())
    },
    
    preprocessWolfData = function() {
      library(gdata)
      
      message("Reading raw data...")
      gps.raw <- read.xls(file.path(context$rawDataDirectory, "GPS_Finland_GPS_Finland_RKTL_Wolf_ec_import_tracks.xlsx"), na.strings="")  
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
      gps.raw <- read.xls(file.path(context$rawDataDirectory, "GPS_Finland_GPS_Finland_RKTL_Lynx_ec_import_tracks.xls"))  
      gps <- data.frame(id=gps.raw$UnitID,
                        date=as.POSIXct(paste(gps.raw$Date, gps.raw$Time), format="%m/%d/%Y %I:%M:%S %p"),
                        x=gps.raw$X, y=gps.raw$Y)
      
      return(list(gps=gps, fromProj="+init=epsg:2393"))
    },
    
    preprocessReindeerData = function() {
      library(gdata)
      
      message("Reading raw data...")
      gps.raw <- read.xls(file.path(context$rawDataDirectory, "GPS_Finland_GPS_Finland_RKTL_ForestReindeer_ec_import_tracks.xlsx"))  
      gps <- data.frame(id=gps.raw$UnitID,
                        date=as.POSIXct(gps.raw$Date*60*60*24-60*60*2*24-60*60*.3, origin="1900-01-01"),
                        x=gps.raw$X, y=gps.raw$Y)
      
      return(list(gps=gps, fromProj="+init=epsg:4326"))
    }
  )
)
