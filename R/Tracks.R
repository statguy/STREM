library(adehabitatLT)
setOldClass("ltraj")

Tracks <- setRefClass(
  Class = "Tracks",
  fields = list(
    study = "Study",
    tracks = "ANY",
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
    
    toGGDF = function(response) {
      library(ggplot2)
      tracksDF <- ggplot2::fortify(getSpatialLines())
      if (!missing(response)) tracksDF$response <- study$getPrettyResponse(response)
      return(tracksDF)
    },
    
    getTracksDirectory = function() return(study$context$processedDataDirectory),
    
    getTracksFileName = function() {
      return(study$context$getFileName(getTracksDirectory(), name="Tracks", response=study$response, region=study$studyArea$region))
    },
    
    saveTracks = function(fileName=getTracksFileName()) {
      stop("Override saveData() method.")
    },
    
    loadTracks = function(fileName=getTracksFileName(), addColumns=TRUE, findTruePopulationSize=TRUE) {
      library(adehabitatLT)
      library(data.table)
      message("Loading tracks from ", fileName, "...")
      load(fileName, envir=as.environment(.self))
      
      if (is.data.frame(tracks)) {
        if (any(!c("year","yday","burst") %in% names(tracks))) {
          message("Breaking down dates...")
          tracks <<- data.frame(tracks, breakDownDate(tracks$date))
          tracks$burst <<- with(tracks, paste(id,year))
        }
        if (addColumns) {
          if (any(!c("dt","dist") %in% names(tracks))) {
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
    
    #getSpatialLines = function(variables=.(id,year)) {
    getSpatialLines = function(variables=.(id,year,burst)) {
      library(sp)
      
      message("Converting tracks to SP object...")
      
      x <- if (inherits(tracks, "ltraj")) {
        x <- ld(tracks)
        x$year <- as.POSIXlt(x$date)$year + 1900
        x
      } else tracks

      if (is.data.frame(x)) {
        library(plyr)
        lines <- dlply(x, variables, function(x, variables) {
          return(Lines(Line(x[,c("x","y")]), ID=paste(x[1,variables], collapse=" ")))
        }, variables=as.character(variables))
        id <- ldply(lines, function(x) data.frame(ID=x@ID))
        tracksSP <- SpatialLinesDataFrame(SpatialLines(lines, proj4string=study$studyArea$proj4string), data=id, match.ID=FALSE)
      }

      return(tracksSP)
    },
    
    plotTracks = function(surveyRoutes, intersects=FALSE, habitat=FALSE) {
      library(sp)
      #library(adehabitatLT)
      
      if (habitat) plot(study$studyArea$habitat)
      
      op <- par(mar=rep(0, 4))
      tracksSP <- getSpatialLines(variables=.(id, year))
      plot(tracksSP, col=1:16, add=habitat)
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
      #distance <<- mean(distances, na.rm=T)
      return(invisible(distances))
    },
    
    getDistances = function(fromSample=TRUE) {
      #library(plyr)
      library(data.table)
      message("Finding distances...")

      tracksDT <- if (inherits(tracks, "ltraj")) data.table(addDtDist(ld(tracks))) else {
        if (fromSample) data.table(tracks[1:10000,]) # TODO: taking random sample would be better, but it has to be taken without breaking the days
        else data.table(tracks)
      }
      
      .internal.getDistances <- function(dt,dist) {
        add <- if (is.na(last(dt))) mean(dt, na.rm=T) else 0
        if (is.na(add)) return(as.numeric(NA))
        s <- (sum(dt, na.rm=T) + add) / 3600
        if (s < 23 | s > 25) return(as.numeric(NA))
        noNAIndex <- !(is.na(dist) | is.na(dt))
        if (length(noNAIndex) == 0) return(as.numeric(NA))
        distances <- sum(dist[noNAIndex]) / sum(dt[noNAIndex]) * 24 * 3600
        return(distances)
      }
      
      tracksDT[, distance:=.internal.getDistances(dt,dist), by=c("year","id","yday")]
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
      #distance <<- mean(distances, na.rm=T)
      #distances <<- rep(distance, times=length(surveyRoutes$surveyRoutes))
      return(invisible(distances))
    },

    getMeanDistance = function() {
      return(mean(getDistances(), na.rm=TRUE))
    },
       
    .internal.thin = function(by) {
      library(plyr)
      
      tracksDF <- if (is.data.frame(tracks)) tracks else ld(tracks)
      #tracksDF <- data.frame(tracksDF, breakDownDate(tracksDF$date))
      oldDt <- mean(tracksDF$dt, na.rm=T)
      
      #tracksDFThinned <- ddply(tracksDF, .variables=.(id, burst, year, yday), .fun=function(x, by) {
      tracksDFThinned <- ddply(tracksDF, .variables=.(id, year, yday), .fun=function(x, by) {
        n <- nrow(x)
        if (n == 1 | n < 2 * by) return(NULL)
        retainIndex <- seq(1, n, by=by)
        return(x[retainIndex,])
      }, by=by, .parallel=TRUE)
      
      if (nrow(tracksDFThinned) == 0) return(NULL)
      
      tracksDFThinned <- tracksDFThinned[,c("x","y","id","burst","date","dt","dist","herdSize","year","month","day","yday","dx","dy")]
      tracksDFThinned <- addDtDist(tracksDFThinned)
      
      newDt <- mean(tracksDFThinned$dt, na.rm=T)
      message("Thinned movements from dt = ", oldDt, " to dt = ", newDt, ", ratio = ", newDt/oldDt)
      return(tracksDFThinned)
    },
    
    getHabitatPreferences = function(habitatWeightsTemplate, nSamples=30, save=FALSE) {
      stop("Unimplemented method.")
    }
  )
)

SimulatedTracks <- setRefClass(
  Class = "SimulatedTracks",
  contains = "Tracks",
  fields = list(
    iteration = "integer",
    truePopulationSize = "data.frame"
  ),
  methods = list(
    initialize = function(xy, id, date, dt, dist, burst, year, yday, herdSize, preprocessData=FALSE, ...) {
      if (!missing(xy) && !missing(id) && !missing(date)) setTracks(xy=xy, id=id, date=date, dt=dt, dist=dist, burst=burst, year=year, yday=yday, herdSize=herdSize)
      callSuper(preprocessData=preprocessData, ...)
    },
    
    setTracks = function(xy, id, date, dt, dist, burst, year, yday, herdSize) {
      #tracks <<- data.frame(xy, id=id, burst=burst, date=date, year=year, yday=yday, dt=dt, dist=dist, herdSize=herdSize)
      tracks <<- data.frame(xy, id=id, burst=burst, date=date, dt=dt, dist=dist, herdSize=herdSize)
      tracks <<- cbind(tracks, breakDownDate(tracks$date))
      tracks <<- addDtDist(tracks)
      return(invisible(.self))
    },
    
    getTracksDirectory = function() return(study$context$scratchDirectory),
    
    getTracksFileIterations = function() {
      if (inherits(study, "undefinedField"))
        stop("Provide study parameter.")
      return(study$context$getIterationIds(dir=getTracksDirectory(), name="Tracks", response=study$response, region=study$studyArea$region, tag="(\\d+)"))
    },
    
    getTracksFileName = function() {
      if (inherits(study, "undefinedField") || length(iteration) == 0)
        stop("Provide study and iteration parameters.")
      return(study$context$getLongFileName(dir=getTracksDirectory(), name="Tracks", response=study$response, region=study$studyArea$region, tag=iteration))
    },
    
    loadTracks = function(fileName=getTracksFileName(), addColumns=TRUE, findTruePopulationSize=TRUE) {
      callSuper(fileName=fileName, addColumns=addColumns)
      if (findTruePopulationSize && (inherits(truePopulationSize, "uninitializedField") || nrow(truePopulationSize)) == 0)
        setTruePopulationSize()
      return(invisible(.self))
    },
    
    saveTracks = function(fileName=getTracksFileName()) {
      message("Saving tracks to ", fileName)
      if (thinId != 1)
        stop("Saving thinned tracks unsupported.")
      save(tracks, iteration, thinId, truePopulationSize, file=fileName)
      return(invisible(.self))
    },
    
    randomizeObservationDayTracks = function(days=1) {
      library(plyr)
      message("Randomizing observation days = ", days, " and filtering tracks...")

      observationTracksDF <- ddply(tracks, .(year, id), function(x, days) {
        days <- days + 1 # Cos we want also the last vector of the day
        randomDays <- x$yday
        if (max(randomDays) + 1 - days < 0)
          stop("Too many days to randomize, max days available = ", max(randomDays) + 1   -1) # So max avail days is -1
        
        randomDays <- randomDays[randomDays <= max(randomDays) + 1 - days]
        randomDay <- base::sample(randomDays, 1)
        
        # From the last day, remove everything except the first relocation
        y <- subset(x, yday %in% randomDay:(randomDay + days - 1))
        lastDay <- max(y$yday)
        y <- rbind(y[y$yday != lastDay,,drop=FALSE], y[y$yday == lastDay,,drop=FALSE][1,,drop=FALSE])
        return(y)
      }, days=days)
      
      observationTracks <- copy(shallow=TRUE) # Avoid possible recursion
      observationTracks$tracks <- observationTracksDF

      return(observationTracks)
    },
    
    getObservationDaysTracks = function(ydays) {
      library(plyr)
      
      observationTracksDF <- ddply(tracks, .(year), function(x, ydays) {
        return(subset(x, yday %in% ydays))
      }, ydays=ydays)
      
      observationTracks <- copy(shallow=TRUE)
      observationTracks$tracks <- observationTracksDF
      return(observationTracks)
    },
    
    getObservationTracksFileName = function() {
      if (inherits(study, "undefinedField") | length(iteration) == 0) 
        stop("Provide study and iteration parameters.")
      return(study$context$getLongFileName(dir=getTracksDirectory(), name="ObservationTracks",
                                           response=study$response, region=study$studyArea$region, tag=iteration))
    },
    
    countIntersections = function(surveyRoutes, days=1, save=TRUE) {
      observationTracks <- randomizeObservationDayTracks(days=days)
      #surveyRoutes <- study$loadSurveyRoutes()
      intersections <- SimulatedIntersections$new(study=study, iteration=iteration)
      intersections$findIntersections(observationTracks, surveyRoutes, dimension=1)
      if (save) {
        intersections$saveIntersections()
        observationTracks$saveTracks(getObservationTracksFileName())
      }
      return(invisible(list(intersections=intersections, tracks=observationTracks, surveyRoutes=surveyRoutes)))
    },
    
    thin = function(by, thinId) {
      message("Thinning iteration = ", iteration, ", thin = ", .self$thinId, " to ", thinId)
      thinnedTracksDF <- .internal.thin(by=by)
      if (is.null(thinnedTracksDF)) return(NULL)
      thinnedTracks <- SimulatedTracks$new(study=study, tracks=thinnedTracksDF, iteration=iteration, thinId=thinId)#thinId=as.integer(thinId + 1))
      return(thinnedTracks)
    },
    
    setTruePopulationSize = function() {
      library(plyr)
      message("Finding true population size...")
      truePopulationSize <<- ddply(tracks, .(year), function(x) {
        n <- daply(x, .(id), function(y) y$herdSize[1])
        return(sum(n))
      })
      
      #return(data.frame(Observed=length(unique(x$id)) * x$herdSize)))
      names(truePopulationSize) <<- c("Year","Observed")
      return(invisible(.self))
    },
    
    getTruePopulationSize = function() {
      if (nrow(truePopulationSize) == 0)
        setTruePopulationSize()
      return(truePopulationSize)
    },
    
    sample = function(nSamples) {
      library(plyr)
      tracks <<- ddply(tracks, .(year), function(x, nSamples) {
        nIds <- length(unique(tracks$id))
        if (nIds < nSamples) nSamples <- nIds
        sampledIds <- base::sample(1:nIds, nSamples)
        return(subset(x, id %in% sampledIds))
      }, nSamples=nSamples)
    },
    
    getHabitatPreferences = function(habitatWeightsTemplate, nSamples=30, save=FALSE) {
      movementIntervals <- ConstantMovementSampleIntervals$new(study=study) # TODO: UNTESTED CODE
      movementIntervals$findSampleIntervals(tracks=.self)
      habitatPreferences <- SimulationHabitatSelection$new(study=study, iteration=iteration)
      habitatPreferences$getHabitatPreferences(intervals=movementIntervals, habitatWeightsTemplate=habitatWeightsTemplate, nSamples=nSamples, save=save)
      return(habitatPreferences)
    }
  )
)

FinlandWTCTracks <- setRefClass(
  Class = "FinlandWTCTracks",
  contains = "Tracks",
  fields = list(
    metadata = "data.frame"
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
      message("Number of individuals in raw data = ", length(unique(gps$id)))
      
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
      message("Number of individuals in raw data = ", length(unique(gps$id)))
      
      return(list(gps=gps, fromProj="+init=epsg:2393"))
    },
    
    preprocessReindeerData = function() {
      library(gdata)
      
      gps.raw <- read.xls(file.path(study$context$rawDataDirectory, "GPS_Finland_GPS_Finland_RKTL_ForestReindeer_ec_import_tracks.xlsx"))
      gps <- data.frame(id=gps.raw$UnitID,
                        date=as.POSIXct(gps.raw$Date*60*60*24-60*60*2*24-60*60*.3, origin="1900-01-01"),
                        x=gps.raw$X, y=gps.raw$Y)
      message("Number of individuals in raw data = ", length(unique(gps$id)))
      
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
    },
    
    preprocessWolfMetadata = function() {
      metadata.raw <- read.csv(file.path(study$context$rawDataDirectory, "wolf-metadata.csv"))
      gps.metadata <- data.frame(id=metadata.raw$nimi,
                                 sex=ifelse(metadata.raw$sp=="uros", "male", "female"),
                                 born=metadata.raw$"syntym\u00e4",
                                 dispersing=metadata.raw$dispersoiva.nuori==1,
                                 homerange=metadata.raw$reviiri==1,
                                 alpha=metadata.raw$alfa==1)
      return(gps.metadata)
    },

    preprocessLynxMetadata = function() {
      gps.metadata <- data.frame(id=c(),
                                 sex=c(),
                                 born=c())
      return(gps.metadata)
    },
    
    preprocessReindeerMetadata = function() {
      gps.metadata <- data.frame(id=c(),
                                 sex=c(),
                                 born=c())
      return(gps.metadata)
    },
    
    getMetadataFileName = function() {
      return(study$context$getFileName(dir=study$context$processedDataDirectory, name="TracksMetadata", response=study$response, region=study$studyArea$region))
    },
    
    saveMetadata = function(fileName=getMetadataFileName()) {
      gps.metadata <- if (study$response == "canis.lupus") preprocessWolfMetadata()
      else if (study$response == "lynx.lynx") preprocessLynxMetadata()
      else if (study$response == "rangifer.tarandus.fennicus") preprocessReindeerMetadata()
      metadata <<- gps.metadata
      save(metadata, file=fileName)
      return(invisible(.self))
    },
    
    loadMetadata = function(fileName=getMetadataFileName()) {
      load(file=fileName, envir=as.environment(.self))
      return(invisible(.self))
    },
    
    #getSpatialLines = function(variables=.(id,year)) {
    getSpatialLines = function(variables=.(id,year,burst)) {
      library(plyr)
      loadMetadata()
      tracksSP <- callSuper(variables=variables)
      if (nrow(metadata) > 0)
        tracksSP@data <- plyr::join(tracksSP@data, metadata, by="id")
      return(tracksSP)
    },
    
    getHabitatPreferences = function(habitatWeightsTemplate, nSamples=30, save=FALSE) {
      movementIntervals <- MaxApproximateConstantMovementSampleIntervals$new(study=study)
      movementIntervals$findSampleIntervals(tracks=.self$tracks)
      habitatPreferences <- WTCHabitatSelection$new(study=study)
      habitatPreferences$getHabitatPreferences(intervals=movementIntervals, habitatWeightsTemplate=habitatWeightsTemplate, nSamples=nSamples, save=save)
      return(habitatPreferences)
    }
  )
)
