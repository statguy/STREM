library(sp)

SurveyRoutes <- setRefClass(
  Class = "SurveyRoutes",
  fields = list(
    study = "Study",
    surveyRoutes = "SpatialLines",
    centroids = "SpatialPoints"
  ),
  methods = list(    
    plotSurveyRoutes = function() {
      plot(study$studyArea$boundary)
      plot(surveyRoutes, col="blue", add=T)
    }
  )
)

FinlandRandomWTCSurveyRoutes <- setRefClass(
  Class = "FinlandRandomWTCSurveyRoutes",
  contains = "SurveyRoutes",
  fields = list(
  ),
  methods = list(
    initialize = function(...) {
      callSuper(...)
      return(.self)
    },
    
    newInstance = function(nStudyRoutes) {
      surveyRoutes <<- getRandomSurveyRoutes(nStudyRoutes=nStudyRoutes)
      return(.self)
    },
    
    getRandomSurveyRoutes = function(nStudyRoutes) {
      initialPopulation <- RandomInitialPopulation$new(studyArea=study$studyArea)
      centroids <<- initialPopulation$randomize(nStudyRoutes)
      angles <- runif(length(centroids), 0, 2*pi)
      return(getTriangles(centroids, angles, 4000))
    }
  )  
)

FinlandWTCSurveyRoutes <- setRefClass(
  Class = "FinlandSurveyRoutes",
  contains = "SurveyRoutes",
  fields = list(
    intersections = "Intersections"  
  ),
  methods = list(
    initialize = function(...) {
      callSuper(...)
      return(.self)
    },
    
    newInstance = function() {
      intersections <<- FinlandWTCIntersections$new(study=study)$newInstance()
      surveyRoutes <<- getWTCSurveyRoutes()
      return(.self)
    },
        
    getWTCSurveyRoutes = function() {
      centroids <<- intersections$getSampleLocations()
      angles <- runif(length(centroids), 0, 2*pi) # Angles are unknown, randomize
      return(getTriangles(centroids, angles, 4000))
    }
  )
)
