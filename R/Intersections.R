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
    
    getIntersectionsFileName = function() {
      if (inherits(study, "undefinedField"))
        stop("Provide study parameters.")
      return(study$context$getFileName(dir=study$context$resultDataDirectory, name="Intersections", response=study$response, region=study$studyArea$region))
    },
    
    saveIntersections = function() {
      stop("Override saveIntersections() method.")
    },
        
    loadIntersections = function() {
      load(getIntersectionsFileName(), envir=as.environment(.self))
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
      #if (length(tracks$distance) == 0)
      #  stop("Did you forgot to run determineDistances() for tracks first?")
      
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
        x <- data.frame(surveyRoute=rownames(intersectionsMatrix[,bursts,drop=F]),
                        x=centroids[,1], y=centroids[,2],
                        year=year,
                        response=study$response,
                        intersections=rowSums(intersectionsMatrix[,bursts,drop=F]),
                        duration=duration,
                        length=surveyRoutes$lengths,
                        distance=1) # Correct for distance after estimation
        data <- rbind(data, x)
      }
      
      intersections <<- SpatialPointsDataFrame(coords=data[,c("x","y")],
                                               data=data[,!names(data) %in% c("x","y")],
                                               proj4string=surveyRoutes$surveyRoutes@proj4string,
                                               match.ID=FALSE)
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
  fields = list(
    covariates = "data.frame"
  ),
  methods = list(
    newInstance = function(...) {
      callSuper(...)
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
      intersections <<- SpatialPointsDataFrame(cbind(wtc$x, wtc$y), data=wtc, proj4string=CRS("+init=epsg:2393"))

      save(intersections, file=getIntersectionsFileName())
    },
    
    getPopulationDensityCovariates = function() {
      library(plyr)
      library(ST)
      
      x <- as.data.frame(intersections)
      populationDensity <- ddply(x, .(year), function(x) {
        year <- x$year[1]
        message("Processing year ", year, "...")
        xy <- SpatialPoints(cbind(x$x, x$y), proj4string=study$studyArea$proj4string)
        populationDensity <- getStatFiPopulationDensity(xy=xy, year=year, aggregationFactor=4) # Grid is 1 x 1 km^2 => aggregate to the scale of the survey routes
        return(data.frame(id=x$id, year=year, populationDensity=populationDensity))
      }, .parallel=FALSE) # Parallel and file caching may cause problems
      
      return(populationDensity)
    },
    
    getWeather = function(year, apiKey) {
      library(ST)
      library(plyr)
      
      query <- "fmi::observations::weather::daily::multipointcoverage"
      startTime <- paste(year, "-01-01T00:00:00Z", sep="")
      endTime <- paste(year, "-03-31T23:59:59Z", sep="")
      parameters <- list(startTime=startTime, endTime=endTime, bbox="19.0900,59.3000,31.5900,70.130", crs="EPSG::4326")
      weatherStationFile <- queryFMIData(apiKey=apiKey, queryStored=query, parameters=parameters)
      weather <- parseFMIWeatherStationMultipointCoverage(weatherStationFile, study$studyArea$proj4string)
      weather <- as.data.frame(weather)
      weather$year <- as.POSIXlt(weather$date)$year + 1900
      weather$month <- as.POSIXlt(weather$date)$mon + 1
      weather$day <- as.POSIXlt(weather$date)$mday
      
      weatherAggregated <- ddply(weather, .(x, y, month), function(x) {
        data.frame(x=x$x[1], y=x$y[1], year=x$year[1], month=x$month[1], rrday=mean(x$rrday, na.rm=T), tday=mean(x$tday, na.rm=T), snow=mean(x$snow, na.rm=T), tmin=mean(x$tmin, na.rm=T), tmax=mean(x$tmax, na.rm=T))
      })
                 
      return(weather)
    },
    
    getWeatherCovariates = function(fmiApiKey) {
      library(plyr)
      
      getWeatherValues = function(weather, variable, x, templateRaster, transform=identity, inverseTransform=identity) {
        weatherRaster <- multiRasterInterpolate(xyz=weather[,c("x","y",variable,"month")],
          variables=.(month),
          templateRaster=templateRaster,
          transform=transform,
          inverseTransform=inverseTransform)        
        
        x$month <- as.POSIXlt(x$date)$mon + 1
        weatherValues <- ddply(x, .(month), function(x, weatherRaster) {
          if (x$month[1] < 1 | x$month[1] > 3) {
            warning("Invalid month ", x$month[1], ", will use February.")
            x$month <- 2
          }
          weatherValues <- extract(weatherRaster[[x$month[1]]], x[,c("x","y")])
          return(data.frame(weatherValues=weatherValues))
        }, weatherRaster=weatherRaster)
        
        return(weatherValues$weatherValues)
      }
      
      templateRaster <- raster(extent(study$studyArea$habitat), nrows=1300, ncols=800)
      x <- as.data.frame(intersections)
      
      weatherCovariates <- ddply(x, .(year), function(x, templateRaster, fmiApiKey) {
        year <- x$year[1]
        message("Processing year ", year, "...")
        
        weather <- getWeather(year=year, apiKey=fmiApiKey)
        snow <- getWeatherValues(weather=weather, variable="snow", x=x, templateRaster=templateRaster, transform=sqrt, inverseTransform=function(x) x^2)
        
        return(data.frame(id=x$id, year=year, snow=snow))
      }, templateRaster=templateRaster, fmiApiKey=fmiApiKey, .parallel=TRUE)
      
      return(weatherCovariates)
    },
    
    getCovariatesFileName = function() {
      if (inherits(study, "undefinedField"))
        stop("Provide study parameters.")
      return(study$context$getFileName(dir=study$context$resultDataDirectory, name="IntersectionsCovariates", response=study$response, region=study$studyArea$region))
    },
    
    saveCovariates = function(fmiApiKey) {
      if (missing(fmiApiKey))
        stop("Provide fmiApiKey argument.")      
      weatherCovariates <- getWeatherCovariates(fmiApiKey=fmiApiKey)
      populationDensityCovariates <- getPopulationDensityCovariates()
      covariates <<- merge(populationDensityCovariates, weatherCovariates, sort=FALSE)
      save(covariates, file=getCovariatesFileName())
      return(invisible(.self))
    },
    
    loadCovariates = function() {
      load(getCovariatesFileName(), envir=as.environment(.self))
      return(invisible(.self))
    }
  )
)
