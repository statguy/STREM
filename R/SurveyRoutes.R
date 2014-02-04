library(sp)

SurveyRoutes <- setRefClass(
  Class = "SurveyRoutes",
  fields = list(
    study = "Study",
    surveyRoutes = "SpatialLines",
    centroids = "SpatialPoints",
    lengths = "numeric"
  ),
  methods = list(    
    plotSurveyRoutes = function() {
      plot(study$studyArea$boundary)
      plot(surveyRoutes, col="blue", add=T)
    },
    
    getLengths = function() {
      library(rgeos)
      library(plyr)
      
      message("Finding survey route lengths...")
      
      contains <- gContains(study$studyArea$boundary, surveyRoutes, byid=T)
      crosses <- gCrosses(study$studyArea$boundary, surveyRoutes, byid=T)
      intersection <- gIntersection(study$studyArea$boundary,
                                    surveyRoutes[c(which(contains), which(crosses))],
                                    byid=T,
                                    id=names(surveyRoutes[c(which(contains), which(crosses))]))
      
      lengths <<- laply(seq(along=intersection), function(i, intersection) {
        coords <- intersection[i]@lines[[1]]@Lines[[1]]@coords
        n <- nrow(coords)
        return(sum(euclidean(coords[1:(n-1),,drop=F], coords[2:n,,drop=F])))
      }, intersection=intersection)
      names(lengths) <<- names(intersection)
      
      outsideNames <- which(!contains & !crosses)
      outsideLengths <- rep(0, length(outsideNames))
      names(outsideLengths) <- outsideNames
      lengths <<- c(lengths, outsideLengths)
      lengths <<- lengths[order(as.numeric(names(lengths)))]
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
      getLengths()
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

## TODO: untested
FinlandWTCSurveyRoutes <- setRefClass(
  Class = "FinlandSurveyRoutes",
  contains = "SurveyRoutes",
  fields = list(
  ),
  methods = list(
    initialize = function(...) {
      callSuper(...)
      return(.self)
    },
    
    newInstance = function() {
      getWTCSurveyRoutes()
      getLengths()
      return(.self)
    },
        
    getWTCSurveyRoutes = function(intersections) {
      intersections <- FinlandWTCIntersections$new(study=study)$newInstance()
      centroids <<- unique(intersections$getSampleLocations())
      angles <- runif(length(centroids), 0, 2 * pi) # Angles are unknown, randomize. TODO: find angles conforming the landscape
      surveyRoutes <<- getTriangles(centroids, angles, 4000)
    }
  )
)
