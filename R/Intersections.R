library(sp)

Intersections <- setRefClass(
  Class = "Intersections",
  contains = "CovariatesContainer",
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
    setDistance = function(distance) {
      intersections$distance <<- distance
      return(invisible(.self))
    },
    
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
      message("Loading intersections from ", fileName, "...")
      load(fileName, envir=as.environment(.self))
      return(invisible(.self))
    },
    
    getSurveyLocations = function() {
      return(SpatialPoints(unique(coordinates(intersections)), proj4string=intersections@proj4string))
    },
    
    getExtent = function(margin=0) {
      library(sp)
      
      coords <- coordinates(intersections)
      index <- intersections$intersections > 0
      x <- range(coords[index,1]) + c(-1,1)*margin
      y <- range(coords[index,2]) + c(-1,1)*margin
      extent <- getPolygonRectangle(x, y, intersections@proj4string)
      
      return(extent)
    },
    
    plotExtent = function(margin=0) {
      extent <- getExtent(margin=margin)
      index <- intersections$intersections > 0
      plot(intersections)
      plot(intersections[index,], col="red", add=T)
      plot(extent, border="blue", add=T)
      return(invisible(.self))
    },
    
    delete = function(extent, permanently=FALSE) {
      index <- over(intersections, extent)
      if (permanently) intersections <<- intersections[!is.na(index),]
      else intersections$intersections[is.na(index)] <<- 0
      return(invisible(.self))
    },
    
    ###
    saveCovariates = function(fileName=getCovariatesFileName()) {
      save(intersections, covariateNames, file=fileName)
      return(invisible(.self))
    },
    
    setCovariatesId = function(tag="") {
      covariatesId <<- paste0("Intersections-", tag)
      return(invisible(.self))
    },
    
    getSampleLocations = function(...) {
      return(intersections[,"date",drop=FALSE])
      return(xyt)
    },
    
    associateCovariates = function(...) {
      covariates <- cbind(...)
      intersections@data <<- cbind(intersections@data, covariates)
      return(invisible(.self))
    },
    
    getCovariates = function() return(intersections@data)
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
      surveyRoutesSP <- surveyRoutes$surveyRoutes
      message("Counting intersections...")
      findIntersectionsMatrix(tracksSP, surveyRoutesSP, ...)
      nSurveyRoutes <- length(surveyRoutes$surveyRoutes)
      nTracks <- length(tracksSP)
      nYears <- length(unique(tracks$tracks$year))
      message("Found ", sum(intersectionsMatrix) / nTracks / nYears, " intersections per track per year.")
      message("Found ", sum(intersectionsMatrix) / nSurveyRoutes / nYears, " intersections per survey route per year.")
      aggregateIntersectionsMatrix(tracks, surveyRoutes)
    },
    
    # TODO: Support for WTC survey routes for each year
    findIntersectionsMatrix = function(tracksSP, surveyRoutesSP, dimension=1) {
      library(sp)
      
      nSurveyRoutes <- length(surveyRoutesSP)
      nTracks <- length(tracksSP)
      #intersectionsMatrix <<- matrix(0, nrow=nSurveyRoutes, ncol=nTracks)
      
      if (dimension == 1) {
        # Assumes that the survey routes are the same every year
        countIntersections <- function(i, surveyRoutesSP, tracksSP) {
          library(plyr)
          #library(rgeos)

          nSurveyRoutes <- length(surveyRoutesSP)
          nTracks <- length(tracksSP)
          cat("Finding intersections for survey route ", i, " / ", nSurveyRoutes, " for ", nTracks, " tracks... ")
                   
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
          
          surveyRouteCoords <- coordinates(surveyRoutesSP@lines[[i]])[[1]]
          
          x <- integer(nTracks)
          for (j in 1:nTracks) {
            trackCoords <- coordinates(tracksSP@lines[[j]])[[1]]
            x[j] <- internal_countIntersections(trackCoords, surveyRouteCoords)
          }
          
          cat("found = ", sum(x), "\n")

          return(x)
        }
        
        intersectionsMatrix <<- laply(1:nSurveyRoutes, countIntersections, surveyRoutesSP=surveyRoutesSP, tracksSP=tracksSP, .drop=FALSE, .parallel=TRUE, .inform=FALSE)
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
      
      #rownames(intersectionsMatrix) <<- names(surveyRoutesSP)
      #colnames(intersectionsMatrix) <<- sapply(tracksSP@lines, function(x) x@ID)

      rownames(intersectionsMatrix) <<- getSPID(surveyRoutesSP)
      colnames(intersectionsMatrix) <<- getSPID(tracksSP)
    },
    
    aggregateIntersectionsMatrix = function(tracks, surveyRoutes) {
      #if (length(tracks$distance) == 0)
      #  stop("Did you forgot to run getDistances() for tracks first?")
      
      distances <- tracks$getDistances()
      distance <- mean(distances, na.rm=T)
      
      tracksObj <- tracks$getTracks()
      tracksDF <- if (inherits(tracksObj, "ltraj")) ld(tracksObj) else tracksObj
      
      burstYear <- ddply(tracksDF, .(burst, id, year), function(x) {
        date <- as.POSIXlt(x$date)
        if (length(unique(x$herdSize)) != 1) stop("Herd size must be constant within a track.")
        y <- data.frame(burst=paste(x[1,c("burst","id","year")], collapse=" "),
                        #burst=x$burst[1],
                        year=date$year[1] + 1900,
                        duration=round(as.numeric(difftime(max(x$date), min(x$date), units="days"))), # Here it is assumed that the observation period is continuous
                        herdSize=x$herdSize[1]) # Herd size does not change within track
        return(y)
      })
      
      data <- data.frame()
      for (year in sort(unique(burstYear$year))) {
        yearToBurstsIndex <- burstYear$year == year
        bursts <- as.character(burstYear[yearToBurstsIndex,]$burst)
        duration <- burstYear[yearToBurstsIndex,]$duration[1]
        herdSize <- burstYear[yearToBurstsIndex,]$herdSize
        centroids <- coordinates(surveyRoutes$centroids)
        x <- data.frame(surveyRoute=as.integer(rownames(intersectionsMatrix[,bursts,drop=F])),
                        x=centroids[,1], y=centroids[,2],
                        year=year,
                        response=study$response,
                        intersections=rowSums(intersectionsMatrix[,bursts,drop=F] * herdSize),
                        duration=duration,
                        length=surveyRoutes$lengths,
                        distance=distance)
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
    },
    
    getIntersectionFileIterations = function() {
      if (inherits(study, "undefinedField"))
        stop("Provide study parameter.")
      return(study$context$getIterationIds(dir=study$context$resultDataDirectory, name="Intersections", response=study$response, region=study$studyArea$region, tag="(\\d+)"))
    }
  )
)

FinlandWTCIntersections <- setRefClass(
  Class = "FinlandWTCIntersections",
  contains = c("Intersections"),
  fields = list(
    maxDuration = "numeric"
  ),
  methods = list(
    initialize = function(maxDuration=Inf, ...) {
      callSuper(...)
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
      
      date <- strptime(wtc$date, "%Y%m%d", tz="EET")
      date <- as.POSIXlt(date)
      index <- date$year == 0
      date$year[index] <- 100
      wtc$date <- as.POSIXct(date)
      wtc <- wtc[date$year >= 80 & date$year <= 150,]
      
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
      
      if (study$response == "rangifer.tarandus.fennicus")
        delete(getPolygonRectangle(c(3.31,3.684)*1e6, c(7.2,6.92)*1e6, intersections@proj4string))
      
      save(intersections, file=getIntersectionsFileName())
      return(invisible(.self))
    }
  )
)

RussiaWTCIntersections <- setRefClass(
  Class = "RussiaWTCIntersections",
  contains = c("Intersections"),
  methods = list(
    initialize = function(...) {
      callSuper(...)
      return(invisible(.self))
    },
    
    saveIntersections = function() {
      library(sp)
      library(dplyr)
      
      message("Processing ", study$response, "...")
      
      districts <- study$studyArea$loadDistricts()
      districts <- data.frame(coordinates(districts), RegionID=districts@data$ADM4_ID, District_Lat=districts@data$NAME_LAT)
      names(districts) <- c("x","y",c("RegionID", "District_Lat"))
      wtc <- read.csv(file.path(study$context$rawDataDirectory, "russia", "Regions_Districts_Latin.csv"), sep="\t", dec=",")
      x <- merge(wtc, districts, by=c("RegionID", "District_Lat"), all.x=TRUE, sort=FALSE)
      species <- switch(study$response, canis.lupus="Волк", lynx.lynx="Рысь", rangifer.tarandus.fennicus="Лесной северный олень")
      
      xy <- x %.% filter(Species_RU == species) %.%
        mutate(intersections = Cnt_trs_Forest + Cnt_trs_Field + Cnt_trs_Bog,
               length = (Length_Forest + Length_Field + Length_Bog) * 1000,
               duration = 1,
               area = districts_Area) %.% # TODO: Area should be on the same scale as rest of the variables
        filter(!is.na(x), !is.na(y), !is.na(length), length > 0, !is.na(intersections))
      coordinates(xy) <- ~x+y
      proj4string(xy) <- study$studyArea$proj4string
      intersections <<- xy
      
      save(intersections, file=getIntersectionsFileName())
    }

  )
)

FinlandRussiaWTCIntersections <- setRefClass(
  Class = "FinlandRussiaWTCIntersections",
  contains = c("Intersections"),
  methods = list(
    initialize = function(...) {
      callSuper(...)
      return(invisible(.self))
    },
    
    saveIntersections = function() {
      library(sp)
      library(rgdal)
      library(plyr)
      
      # TODO: divide finland into smaller regions
      
      finlandStudy <- FinlandWTCStudy$new(context=context, response=study$response)
      russiaStudy <- RussiaWTCStudy$new(context=context, response=study$response)
      finland <- FinlandWTCIntersections$new(study=finlandStudy)$loadIntersections()
      finland$intersections <- spTransform(finland$intersections, russiaStudy$studyArea$proj4string)
      x <- ddply(as.data.frame(finland$intersections), .(year),
                 function(x) { data.frame(year=x$year[1], length=sum(x$length), duration=mean(x$duration), intersections=sum(x$intersections)) })
      x$district <- "Finland"
      x$area <- finlandStudy$studyArea$boundary@polygons[[1]]@area / 1000^2
      finlandBoundary <- spTransform(finlandStudy$studyArea$boundary, russiaStudy$studyArea$proj4string)
      x$x <- coordinates(finlandBoundary)[1]
      x$y <- coordinates(finlandBoundary)[2]
      
      russia <- RussiaWTCIntersections$new(study=russiaStudy)$loadIntersections()
      y <- as.data.frame(russia$intersections)[,c("x","y","length","intersections","duration","area")]
      y$year <- russia$intersections$Year
      y$district <- russia$intersections$NAME_LAT
      
      z <- rbind(x[,colnames(y)],y)
      coordinates(z) <- ~x+y
      proj4string(z) <- russiaStudy$studyArea$proj4string
      
      minx <- xmin(extent(study$studyArea$habitat))
      maxx <- xmax(extent(study$studyArea$habitat))
      miny <- ymin(extent(study$studyArea$habitat))
      maxy <- ymax(extent(study$studyArea$habitat))
      window <- matrix(c(minx,miny, minx,maxy, maxx,maxy, maxx,miny, minx,miny), ncol=2, byrow=T)
      coords <- coordinates(z)
      intersections <<- z[point.in.polygon(coords[,1], coords[,2], window[,1], window[,2])==1,]
      
      #Remove samples outside the study area
      #uniqCoords <- unique(coordinates(intersections))
      #r <- study$studyArea$habitat
      #remove <- logical(nrow(uniqCoords))
      #for (row in 1:nrow(uniqCoords)) {
      #  x <- uniqCoords[row,1]
      #  y <- uniqCoords[row,2]
      #  rcol <- colFromX(r, x)
      #  rrow <- rowFromY(r, y)
      #  if (r[rrow, rcol] %in% c(0,255) | is.na(r[rrow, rcol])) remove[row] <- TRUE  
      #}
      #join(intersections$intersections, cbind(uniqCoords, remove)
            
      save(intersections, file=getIntersectionsFileName())
    }
  )
)
