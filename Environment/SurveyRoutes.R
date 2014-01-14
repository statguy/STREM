library(sp)

rotate <- function(x, y, angle) {
  return(cbind(x * cos(angle) - y * sin(angle), x * sin(angle) + y * cos(angle)))
}

getTriangles <- function(centroids, angles, sideLength) {
  require(plyr)
  
  l2 <- sideLength / 2
  lmid <- l2 * sin(pi / 3)
  xy <- coordinates(centroids)
  c1 <- rotate(-l2, lmid, angles) + xy
  c2 <- rotate(l2, lmid, angles) + xy
  c3 <- rotate(0, -lmid, angles) + xy
  
  triangles <- llply(1:nrow(c1), function(i, c1, c2, c3) {
    return(Lines(list(Polygon(rbind(c1[i,], c2[i,], c3[i,], c1[i,]))), ID=i))
  }, c1=c1, c2=c2, c3=c3)
  triangles <- SpatialLines(triangles, proj4string=centroids@proj4string)
  return(triangles)
}

SurveyRoutes <- setRefClass(
  Class = "SurveyRoutes",
  fields = list(
    surveyRoutes = "SpatialLines"
  ),
  methods = list(
    findIntersections = function(tracks, runParallel=TRUE, cluster, dimension) {
      require(sp)
      require(parallel)
      
      if (is.null(surveyRoutes)) return(NULL)
      
      nSurveyRoutes <- length(surveyRoutes)
      nTracks <- length(tracks)
      intersections <- matrix(0, nrow=nSurveyRoutes, ncol=nTracks)
      rownames(intersections) <- names(surveyRoutes)
      colnames(intersections) <- names(tracks@lines)
      
      if (missing(cluster) | missing(dimension)) {
        require(rgeos)
        require(plyr)
        for (j in 1:nTracks) {
          message("Finding intersections for track ", j, " / ", nTracks, "...")
          intersections[,j] <- laply(1:nSurveyRoutes,
            function(i, surveyRoutes, tracks, j) {
              return(length(gIntersection(surveyRoutes[i], tracks[j,], byid=TRUE)))
            }, surveyRoutes=surveyRoutes, tracks=tracks, j=j, .parallel=runParallel)
        }
      }
      else if (dimension == 1) {
        intersections <- parSapply(cluster$remoteCluster, 1:nTracks, function(j, surveyRoutes, tracks, cluster) {
          nSurveyRoutes <- length(surveyRoutes)
          nTracks <- length(tracks)
          message("Finding intersections for track ", j," / ", nTracks, "...")
          
          require(plyr)
          require(parallel)
          require(doMC)
          registerDoMC(cores=detectCores())
          require(rgeos)
          
          x <- laply(1:nSurveyRoutes,
            function(i, surveyRoutes, tracks, j) {
              return(length(gIntersection(surveyRoutes[i], tracks[j,], byid=TRUE)))
            }, surveyRoutes=surveyRoutes, tracks=tracks, j=j, .parallel=TRUE)
          
          return(x)
        }, surveyRoutes=surveyRoutes, tracks=tracks)
      }
      else if (dimension == 2) {
        intersections <- parSapply(cluster$remoteCluster, 1:nSurveyRoutes, function(i, surveyRoutes, tracks) {
          nSurveyRoutes <- length(surveyRoutes)
          nTracks <- length(tracks)
          message("Finding intersections for survey route ", i, " / ", nSurveyRoutes, "...")
          
          require(plyr)
          require(parallel)
          require(doMC)
          registerDoMC(cores=detectCores())
          require(rgeos)
          
          x <- laply(1:nTracks,
            function(j, surveyRoutes, tracks, i) {
              return(length(gIntersection(surveyRoutes[i], tracks[j,], byid=TRUE)))
            }, surveyRoutes=surveyRoutes, tracks=tracks, i=i, .parallel=TRUE)
          
          return(x)
        }, surveyRoutes=surveyRoutes, tracks=tracks)
        
        intersections <- t(intersections)
      }
      
      message("Found ", sum(intersections) / nTracks, " intersections per animal track.")
      message("Found ", sum(intersections) / nSurveyRoutes, " intersections per survey route.")
      
      return(intersections)
    }
  )
)

FinlandRandomSurveyRoutes <- setRefClass(
  Class = "FinlandRandomSurveyRoutes",
  contains = "SurveyRoutes",
  fields = list(
    context = "Context"
  ),
  methods = list(
    initialize = function(...) {
      callSuper(...)
      return(.self)
    },
    
    newInstance = function(nStudyRoutes) {
      studyArea <- FinlandStudyArea(context=context)$newInstance()
      surveyRoutes <<- getRandomSurveyRoutes(nStudyRoutes=nStudyRoutes, studyArea$boundary)
      return(.self)
    },
    
    getRandomSurveyRoutes = function(nStudyRoutes, boundary) {
      initialPopulation <- RandomInitialPopulation(boundary)
      centroids <- initialPopulation$randomize(nStudyRoutes)
      angles <- runif(length(centroids), 0, 2*pi)
      return(getTriangles(centroids, angles, 4000))
    }
    
  )  
)

FinlandWTCSurveyRoutes <- setRefClass(
  Class = "FinlandSurveyRoutes",
  contains = "SurveyRoutes",
  fields = list(
    context = "Context",
    wtcData = "WTCData"  
  ),
  methods = list(
    initialize = function(...) {
      callSuper(...)
      return(.self)
    },
    
    newInstance = function() {
      wtcData <<- FinlandWTCData(context=context)$newInstance()
      surveyRoutes <<- getWTCSurveyRoutes()
      return(.self)
    },
        
    getWTCSurveyRoutes = function() {
      centroids <- wtcData$getSampleLocations()
      angles <- runif(length(centroids), 0, 2*pi) # Angles are unknown, randomize
      return(getTriangles(centroids, angles, 4000))
    }
    
  )
)
