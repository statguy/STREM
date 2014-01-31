library(sp)

WTCData <- setRefClass(
  Class = "WTCData",
  fields = list(
    context = "Context",
    studyArea = "StudyArea",
    data = "SpatialPointsDataFrame"
  ),
  methods = list(
    initialize = function(preprocessData=FALSE, ...) {
      callSuper(...)
      if (preprocessData) saveData()
      return(.self)
    },
    
    newInstance = function() {
      data <<- loadData()
    },
    
    getDataFileName = function() {
      return(context$getFileName(context$processedDataDirectory, name="WTC", region=studyArea$region))
    },
    
    saveData = function() {
      stop("Override saveData() method.")
    },
    
    getSampleLocations = function() {
      return(SpatialPoints(unique(coordinates(data)), proj4string=data@proj4string))
    },
    
    loadData = function() {
      load(getDataFileName())
      return(wtc)
    }
  )
)

FinlandWTCData <- setRefClass(
  Class = "FinlandWTCData",
  contains = "WTCData",
  methods = list(
    newInstance = function() {
      callSuper()
      studyArea <<- FinlandStudyArea(context=context)$newInstance()
      return(.self)
    },
    
    saveData = function() {
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
      save(wtc, file=getDataFileName())
    }
  )
)
