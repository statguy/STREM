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
