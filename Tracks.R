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
      load(getDataFileName())
      return(tracks)
    },
    
    plotTracks = function() {
      require(sp)
      require(adehabitatLT)
      require(plyr)
      
      tracksDF <- ld(tracks)
      plot(tracksDF$x, tracksDF$y, type="n")
      ddply(tracksDF, .(id), function(tracksDF) {
        lines(tracksDF$x, tracksDF$y, col=tracksDF$id)
      })
      plot(study$studyArea$boundary, add=T)
    }
  )
)

SimulatedTracks <- setRefClass(
  Class = "SimulationTracks",
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
      if (class(tracksDF) == "uninitializedField")
        stop("Provide tracksDF argument.")
      
      xy <- tracksDF[,c("x","y")]
      date <- as.POSIXct(strptime(paste(2000+tracksDF$year, tracksDF$day, tracksDF$hour, "0", "0"), format="%Y %j %H %M %S"))
      id <- tracksDF$agent
      tracks <<- as.ltraj(xy=xy, date=date, id=id, infolocs=NULL)
      iteration <<- tracksDF$iteration[1]      
    },
    
    saveTracks = function() {      
      save(tracks, file=getDataFileName())
    },
    
    getDataFileName = function() {
      if (class(study) == "undefinedField" | length(iteration) == 0)
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
      
      require(sp)
      require(plyr)
      require(rgdal)
      require(adehabitatLT)
      
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
      require(gdata)
      
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
      require(gdata)

      message("Reading raw data...")
      gps.raw <- read.xls(file.path(context$rawDataDirectory, "GPS_Finland_GPS_Finland_RKTL_Lynx_ec_import_tracks.xls"))  
      gps <- data.frame(id=gps.raw$UnitID,
                        date=as.POSIXct(paste(gps.raw$Date, gps.raw$Time), format="%m/%d/%Y %I:%M:%S %p"),
                        x=gps.raw$X, y=gps.raw$Y)
      
      return(list(gps=gps, fromProj="+init=epsg:2393"))
    },
    
    preprocessReindeerData = function() {
      require(gdata)
      
      message("Reading raw data...")
      gps.raw <- read.xls(file.path(context$rawDataDirectory, "GPS_Finland_GPS_Finland_RKTL_ForestReindeer_ec_import_tracks.xlsx"))  
      gps <- data.frame(id=gps.raw$UnitID,
                        date=as.POSIXct(gps.raw$Date*60*60*24-60*60*2*24-60*60*.3, origin="1900-01-01"),
                        x=gps.raw$X, y=gps.raw$Y)
      
      return(list(gps=gps, fromProj="+init=epsg:4326"))
    }
  )
)
