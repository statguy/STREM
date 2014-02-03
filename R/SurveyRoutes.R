library(sp)

SurveyRoutes <- setRefClass(
  Class = "SurveyRoutes",
  fields = list(
    studyArea = "StudyArea",
    surveyRoutes = "SpatialLines"
  ),
  methods = list(
    rotate = function(x, y, angle) {
      return(cbind(x * cos(angle) - y * sin(angle), x * sin(angle) + y * cos(angle)))
    },
    
    plotSurveyRoutes = function() {
      plot(studyArea$boundary)
      plot(surveyRoutes, col="blue", add=T)
    }
      
    #findIntersections = function(tracks, runParallel=TRUE, cluster, dimension) {
    #}
  )
)

TriangleSurveyRoutes <- setRefClass(
  Class = "TriangleSurveyRoutes",
  contains = "SurveyRoutes",
  fields = list(
  ),
  methods = list(
    getTriangles = function(centroids, angles, sideLength) {
      library(plyr)
      
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
    
    #findIntersections = function(tracks, runParallel=TRUE, cluster, dimension) {
    #}
  )
)

FinlandRandomWTCSurveyRoutes <- setRefClass(
  Class = "FinlandRandomWTCSurveyRoutes",
  contains = "TriangleSurveyRoutes",
  fields = list(
    context = "Context"
  ),
  methods = list(
    initialize = function(context=context, ...) {
      if (missing(context))
        stop("Provide context.")
      callSuper(context=context, ...)
      return(.self)
    },
    
    newInstance = function(nStudyRoutes) {
      surveyRoutes <<- getRandomSurveyRoutes(nStudyRoutes=nStudyRoutes)
      return(.self)
    },
    
    getRandomSurveyRoutes = function(nStudyRoutes) {
      initialPopulation <- RandomInitialPopulation$new(studyArea=studyArea)
      centroids <- initialPopulation$randomize(nStudyRoutes)
      angles <- runif(length(centroids), 0, 2*pi)
      return(getTriangles(centroids, angles, 4000))
    }
  )  
)

FinlandWTCSurveyRoutes <- setRefClass(
  Class = "FinlandSurveyRoutes",
  contains = "TriangleSurveyRoutes",
  fields = list(
    context = "Context",
    intersections = "Intersections"  
  ),
  methods = list(
    initialize = function(context, ...) {
      if (missing(context))
        stop("Provide context.")
      callSuper(context=context, ...)
      return(.self)
    },
    
    newInstance = function() {
      intersections <<- FinlandWinterTrackCounts$new(context=context)$newInstance()
      surveyRoutes <<- getWTCSurveyRoutes()
      return(.self)
    },
        
    getWTCSurveyRoutes = function() {
      centroids <- intersections$getSampleLocations()
      angles <- runif(length(centroids), 0, 2*pi) # Angles are unknown, randomize
      return(getTriangles(centroids, angles, 4000))
    }
  )
)
