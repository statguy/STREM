library(adehabitatLT)
setOldClass("ltraj")

Tracks <- setRefClass(
  Class = "Tracks",
  fields = list(
    study = "Study",
    tracks = "ANY",
    distance = "ANY",
    thinId = "integer"
  ),
  methods = list(
    initialize = function(thinId, preprocessData=FALSE, ...) {
      callSuper(...)
      if (missing(thinId)) thinId <<- as.integer(1)
      else thinId <<- thinId
      if (preprocessData) saveTracks()
      return(.self)
    },
    
    getTracksDirectory = function() return(study$context$processedDataDirectory),
    
    getTracksFileName = function() {
      return(study$context$getFileName(getTracksDirectory(), name="Tracks", response=study$response, region=study$studyArea$region))
    },
    
    saveTracks = function(fileName=getTracksFileName()) {
      stop("Override saveData() method.")
    },
    
    loadTracks = function(fileName=getTracksFileName(), addColumns=TRUE) {
      library(adehabitatLT)
      load(fileName, envir=as.environment(.self))
      
      if (is.data.frame(tracks)) {
        if (any(!names(tracks) %in% c("year","yday","burst"))) {
          message("Breaking down dates...")
          tracks <<- data.frame(tracks, breakDownDate(tracks$date))
          tracks$burst <<- with(tracks, paste(id,year))
        }
        if (addColumns) {
          if (any(!names(tracks) %in% c("dt","dist"))) {
            message("Finding sample intervals and distances...")
            tracks <<- addDtDist(tracks)
          }
        }
      }
      
      return(invisible(.self))
    },
    
    getTracks = function() {
      if (nrow(tracks) == 0) loadTracks()
      return(tracks)
    },
        
    getSpatialLines = function(variables=.(burst)) {
      library(sp)
      if (is.data.frame(tracks)) {
        library(plyr)
        lines <- dlply(tracks, variables, function(x, variables) {
          return(Lines(Line(x[,c("x","y")]), ID=paste(x[1,variables], collapse=" ")))
        }, variables=as.character(variables))
        tracksSP <- SpatialLines(lines, proj4string=study$studyArea$proj4string)
      }
      else if (inherits(tracks, "ltraj")) {
        tracksSP <- ltraj2sldf(tracks, byid=TRUE)
        proj4string(tracksSP) <- study$studyArea$proj4string
      }
      return(tracksSP)
    },
    
    plotTracks = function(surveyRoutes, intersects=FALSE) {
      library(sp)
      #library(adehabitatLT)
      
      op <- par(mar=rep(0, 4))
      tracksSP <- getSpatialLines(variables=.(id, year))
      plot(tracksSP, col=1:16)
      plot(study$studyArea$boundary, add=T)
      if (!missing(surveyRoutes)) {
        col <- rep("blue", length(intersects))
        col[intersects] <- "red"
        plot(surveyRoutes$surveyRoutes, col=col, add=T)
      }
      par(op)
      
      #tracksDF <- if (is.data.frame(tracks)) tracks else ld(tracks)
      #plot(tracksDF$x, tracksDF$y, type="n", asp=1)
      #apply(fun=function(x) lines(x$x, x$y, col=x$id), variables=.(burst))
      #plot(study$studyArea$boundary, add=T)
      #if (!missing(surveyRoutes))
      #  plot(surveyRoutes$surveyRoutes, col="blue", add=T)
    },
    
    apply = function(variables=.(id), fun, ..., combine=F) {
      tracksDF <- if (is.data.frame(tracks)) tracks else ld(tracks)
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
    NEW_getDistances = function() {
      library(dplyr)
      message("Finding distances...")
      
      tracksDF <- if (inherits(tracks, "ltraj")) addDtDist(ld(tracks)) else tracks
      
      dayDistance <- function(dt, dist) {
        add <- if (is.na(last(dt))) mean(dt, na.rm=T) else 0
        if (is.na(add)) return(as.numeric(NA))
        s <- (sum(dt, na.rm=T) + add) / 3600
        if (s < 23 | s > 25) return(as.numeric(NA))
        noNAIndex <- !(is.na(dist) | is.na(dt))
        if (length(noNAIndex) == 0) return(as.numeric(NA))
        return(sum(dist[noNAIndex]) / sum(dt[noNAIndex]) * 24 * 3600)        
      }
      
      distances <- tracksDF %.%
        group_by(year, id, yday) %.%
        summarise(distance=dayDistance(dt, dist))
      distances <- distances$distance
      
      if (all(is.na(distances)))
        stop("Unable to determine movement distance.")
      
      distances <- distances[!is.na(distances)]
      if (any(distances > 50000)) {
        warning("Unrealistic distances. Removing those...")
        distances <- distances[distances<50000]
      }
      
      message(sum(!is.na(distances)), " / ", length(distances), " of day movements used to determine mean movement distance = ", mean(distances, na.rm=T), " ± ", sd(distances, na.rm=T), ".")
      distance <<- mean(distances, na.rm=T)
      return(invisible(distances))
    },
    
    getDistances = function() {
      #library(plyr)
      library(data.table)
      message("Finding distances...")
      
      #tracksDF <- if (inherits(tracks, "ltraj")) addDtDist(ld(tracks)) else tracks
      tracksDT <- if (inherits(tracks, "ltraj")) data.table(addDtDist(ld(tracks))) else data.table(tracks)
      
      .internal.getDistances <- function(dt,dist) {
        add <- if (is.na(last(dt))) mean(dt, na.rm=T) else 0
        if (is.na(add)) return(NA)
        s <- (sum(dt, na.rm=T) + add) / 3600
        if (s < 23 | s > 25) return(NA)
        noNAIndex <- !(is.na(dist) | is.na(dt))
        if (length(noNAIndex) == 0) return(NA)
        return(sum(dist[noNAIndex]) / sum(dt[noNAIndex]) * 24 * 3600)        
      }
      
      tracksDT[, distance:=.internal.getDistances(dt,dist), by=c(year,id,yday)]
      distances <- tracksDT$distance
      
      #distances <- ddply(tracksDF, .(burst, year, id, yday), function(x) {
      #  add <- if (is.na(x$dt[nrow(x)])) mean(x$dt, na.rm=T) else 0
      #  if (is.na(add)) return(NA)
      #  s <- (sum(x$dt, na.rm=T) + add) / 3600
      #  if (s < 23 | s > 25) return(NA)
      #  noNAIndex <- !(is.na(x$dist) | is.na(x$dt))
      #  if (length(noNAIndex) == 0) return(NA)
      #  return(sum(x$dist[noNAIndex]) / sum(x$dt[noNAIndex]) * 24 * 3600)
      #}, .parallel=TRUE, .inform=TRUE)$V1
      
      if (all(is.na(distances)))
        stop("Unable to determine movement distance.")
      
      distances <- distances[!is.na(distances)]
      if (any(distances > 50000)) {
        warning("Unrealistic distances. Removing those...")
        distances <- distances[distances<50000]
      }
      
      message(sum(!is.na(distances)), " / ", length(distances), " of day movements used to determine mean movement distance = ", mean(distances, na.rm=T), " ± ", sd(distances, na.rm=T), ".")
      distance <<- mean(distances, na.rm=T)
      #distances <<- rep(distance, times=length(surveyRoutes$surveyRoutes))
      return(invisible(distances))
    },

    getMeanDistance = function() {
      return(mean(getDistances(), na.rm=TRUE))
    },
    
    # Thinning strategies:
    # 
    # xxxxxxxxxxx
    # x x x x x x
    # x   x   x
    # x       x
    #
    # xxxxxxxxxx
    # x x x x x
    # x  x  x  x
    # x   x   x
    # x    x   
    
    .internal.thin = function(by) {
      library(plyr)
      
      tracksDF <- if (is.data.frame(tracks)) tracks else ld(tracks)
      #tracksDF <- data.frame(tracksDF, breakDownDate(tracksDF$date))
      oldDt <- mean(tracksDF$dt, na.rm=T)
      
      #tracksDFThinned <- ddply(tracksDF, .variables=.(id, burst, year, yday), .fun=function(x, by) {
      tracksDFThinned <- ddply(tracksDF, .variables=.(burst), .fun=function(x, by) {
        n <- nrow(x)
        if (n == 1 | n < 2*by) return(NULL)
        retainIndex <- seq(1, n, by=by)
        return(x[retainIndex,])
      }, by=by, .parallel=TRUE)
      
      if (nrow(tracksDFThinned) == 0) return(NULL)
      tracksDFThinned <- tracksDFThinned[,c("x","y","id","date","year","yday","burst","dt","dist")]
      tracksDFThinned <- addDtDist(tracksDFThinned)
      
      newDt <- mean(tracksDFThinned$dt, na.rm=T)
      message("Thinned movements from dt = ", oldDt, " to dt = ", newDt, ", ratio = ", newDt/oldDt)
      return(tracksDFThinned)
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
    iteration = "integer"
  ),
  methods = list(
    initialize = function(xy, id, date, dt, dist, burst, year, yday, preprocessData=FALSE, ...) {
      if (!missing(xy) && !missing(id) && !missing(date)) setTracks(xy=xy, id=id, date=date, dt=dt, dist=dist, burst=burst, year=year, yday=yday)
      callSuper(preprocessData=preprocessData, ...)
    },
    
    setTracks = function(xy, id, date, dt, dist, burst, year, yday) {  
      tracks <<- data.frame(xy, id=id, burst=burst, date=date, year=year, yday=yday, dt=dt, dist=dist)
      return(invisible(.self))
    },
    
    getTracksDirectory = function() return(study$context$scratchDirectory),
    
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
      save(tracks, iteration, thinId, distance, file=fileName)
      return(invisible(.self))
    },
    
    randomizeObservationDayTracks = function() {
      library(plyr)
      message("Randomizing observation day and filtering tracks...")

      observationTracksDF <- ddply(tracks, .(year), function(x) {
        randomDay <- sample(x$yday, 1)
        return(subset(x, yday == randomDay))
      })
      
      observationTracks <- copy()
      observationTracks$tracks <- observationTracksDF
      return(observationTracks)
    },
    
    getObservationDaysTracks = function(ydays) {
      library(plyr)
      
      observationTracksDF <- ddply(tracks, .(year), function(x, ydays) {
        return(subset(x, yday %in% ydays))
      }, ydays=ydays)
      
      observationTracks <- copy()
      observationTracks$tracks <- observationTracksDF
      return(observationTracks)
    },
    
    countIntersections = function() {
      observationTracks <- randomizeObservationDayTracks()
      surveyRoutes <- study$loadSurveyRoutes()
      intersections <- SimulatedIntersections$new(study=study, iteration=iteration)
      intersections$findIntersections(observationTracks, surveyRoutes, dimension=1)
      intersections$saveIntersections()
      return(invisible(.self))
    },
    
    thin = function(by, thinId) {
      message("Thinning iteration = ", iteration, ", thin = ", .self$thinId, " to ", thinId)
      thinnedTracksDF <- .internal.thin(by=by)
      if (is.null(thinnedTracksDF)) return(NULL)
      thinnedTracks <- SimulatedTracks$new(study=study, tracks=thinnedTracksDF, iteration=iteration, thinId=thinId)#thinId=as.integer(thinId + 1))
      return(thinnedTracks)
    },
    
    getTruePopulationSize = function() {
      library(plyr)
      populationSize <- ddply(tracks, .(year), function(x) return(data.frame(Observed=length(unique(x$id)))))
      names(populationSize) <- c("Year","Observed")
      return(populationSize)
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
      
      message("Processing ", study$response, "...")
      
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
      
      tracks <<- gpsLT
      save(tracks, file=getTracksFileName())
      return(invisible(.self))
    },
    
    preprocessWolfData = function() {
      library(gdata)
      
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
      
      gps.raw <- read.xls(file.path(study$context$rawDataDirectory, "GPS_Finland_GPS_Finland_RKTL_Lynx_ec_import_tracks.xls"))  
      gps <- data.frame(id=gps.raw$UnitID,
                        date=as.POSIXct(paste(gps.raw$Date, gps.raw$Time), format="%m/%d/%Y %I:%M:%S %p"),
                        x=gps.raw$X, y=gps.raw$Y)
      
      return(list(gps=gps, fromProj="+init=epsg:2393"))
    },
    
    preprocessReindeerData = function() {
      library(gdata)
      
      gps.raw <- read.xls(file.path(study$context$rawDataDirectory, "GPS_Finland_GPS_Finland_RKTL_ForestReindeer_ec_import_tracks.xlsx"))
      gps <- data.frame(id=gps.raw$UnitID,
                        date=as.POSIXct(gps.raw$Date*60*60*24-60*60*2*24-60*60*.3, origin="1900-01-01"),
                        x=gps.raw$X, y=gps.raw$Y)
      
      return(list(gps=gps, fromProj="+init=epsg:4326"))
    },
    
    getSampleIntervals = function() {
      intervals <- FinlandMovementSampleIntervals$new(study=study)
      intervals$getSampleIntervals(tracks=.self)
      return(intervals)
    },
    
    thin = function(by, thinId) {
      message("Thinning thin = ", thinId)
      thinnedTracksDF <- .internal.thin(by=by)
      if (is.null(thinnedTracksDF)) return(NULL)
      thinnedTracks <- FinlandWTCTracks$new(study=study, tracks=thinnedTracksDF, thinId=thinId)#as.integer(thinId+1))
      return(thinnedTracks)
    }
  )
)
