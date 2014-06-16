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
      return(invisible(.self))
    },

    getSurveyRoutesFileName = function() {
      return(context$getFileName(dir=context$resultDataDirectory, name="SurveyRoutes", response=study$response, region=study$studyArea$region))
    },
    
    loadSurveyRoutes = function(fileName=getSurveyRoutesFileName(), findLengths) {
      load(fileName, env=as.environment(.self))
      return(invisible(.self))
    },

    saveSurveyRoutes = function(fileName=getSurveyRoutesFileName(), findLengths=TRUE) {
      save(surveyRoutes, centroids, lengths, file=fileName)
      if (findLengths) getLengths()
      return(invisible(.self))
    },
    
    cutSurveyRoutes = function() {
      library(rgeos)
      
      contains <- gContains(study$studyArea$boundary, surveyRoutes, byid=T)
      crosses <- gCrosses(study$studyArea$boundary, surveyRoutes, byid=T)
      index <- c(which(contains), which(crosses))
      intersection <- gIntersection(study$studyArea$boundary, surveyRoutes[index], byid=T, id=names(surveyRoutes[index]))
      
      return(intersection)
    },
    
    getLengths = function() {
      library(rgeos)
      library(plyr)
      message("Finding survey route lengths...")
      
      contains <- gContains(study$studyArea$boundary, surveyRoutes, byid=T)
      crosses <- gCrosses(study$studyArea$boundary, surveyRoutes, byid=T)
      index <- c(which(contains), which(crosses))
      intersection <- gIntersection(study$studyArea$boundary, surveyRoutes[index], byid=T, id=names(surveyRoutes[index]))
      
      #lengths <<- laply(seq(along=intersection), function(i, intersection) {
      #  coords <- intersection[i]@lines[[1]]@Lines[[1]]@coords
      #  n <- nrow(coords)
      #  return(sum(euclidean(coords[1:(n-1),,drop=F], coords[2:n,,drop=F])))
      #}, intersection=intersection)
        
      lengths <<- laply(seq(along=intersection), function(i, intersection) {
        length <- laply(seq(along=intersection[i]@lines[[1]]@Lines), function(j, lines) {
          #coords <- intersection[i]@lines[[1]]@Lines[[j]]@coords
          coords <- lines[[j]]@coords
          n <- nrow(coords)
          length <- sum(euclidean(coords[1:(n - 1),, drop=F], coords[2:n,, drop=F]))
          return(length)
        }, lines=intersection[i]@lines[[1]]@Lines)
        return(sum(length))
      }, intersection=intersection)
      
      names(lengths) <<- names(intersection)
      
      outsideNames <- which(!contains & !crosses)
      outsideLengths <- rep(0, length(outsideNames))
      names(outsideLengths) <- outsideNames
      lengths <<- c(lengths, outsideLengths)
      lengths <<- lengths[order(as.numeric(names(lengths)))]
      
      return(invisible(.self))
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
      return(invisible(.self))
    },
    
    randomizeSurveyRoutes = function(nSurveyRoutes, save=FALSE) {
      initialPopulation <- RandomInitialPopulation$new(studyArea=study$studyArea)
      centroids <<- initialPopulation$randomize(nSurveyRoutes)
      angles <- runif(length(centroids), 0, 2*pi)
      surveyRoutes <<- getTriangles(centroids, angles, 4000)
      getLengths()
      if (save) saveSurveyRoutes()
      return(invisible(.self))
    },
    
    getSurveyRoutesFileName = function() {
      return(context$getFileName(dir=context$resultDataDirectory, name="SurveyRoutes", response="Random", region=study$studyArea$region))
    }
  )  
)

FinlandWTCSurveyRoutes <- setRefClass(
  Class = "FinlandSurveyRoutes",
  contains = "SurveyRoutes",
  fields = list(
  ),
  methods = list(
    initialize = function(...) {
      callSuper(...)
      return(invisible(.self))
    },
    
    loadSurveyRoutes = function(fileName=getSurveyRoutesFileName(), findLengths=TRUE) {
      intersections <- FinlandWTCIntersections$new(study=study)$loadIntersections()
      centroids <<- intersections$getSurveyLocations()
      # Angles are unknown, randomize. TODO: Find the most likely angles given the landscape.
      angles <- runif(length(centroids), 0, 2 * pi)
      surveyRoutes <<- getTriangles(centroids, angles, 4000)
      if (findLengths) getLengths()
      return(invisible(.self))
    },
    
    saveSurveyRoutes = function(fileName=getSurveyRoutesFileName()) {
      stop("Unsupported method.")
    }
  )
)

TestSurveyRoutes <- setRefClass(
  Class = "TestSurveyRoutes",
  contains = "SurveyRoutes",
  fields = list(
  ),
  methods = list(
    initialize = function(...) {
      callSuper(...)
      return(invisible(.self))
    },
    
    randomizeSurveyRoutes = function(nSurveyRoutes, save=FALSE) {
      initialPopulation <- RandomInitialPopulation$new(studyArea=study$studyArea)
      centroids <<- initialPopulation$randomize(nSurveyRoutes)
      angles <- runif(length(centroids), 0, 2*pi)
      surveyRoutes <<- getTriangles(centroids, angles, 4000)
      getLengths()
      if (save) saveSurveyRoutes()
      return(invisible(.self))
    },
    
    getSurveyRoutesFileName = function() {
      return(context$getFileName(dir=context$resultDataDirectory, name="SurveyRoutes", response="Random", region=study$studyArea$region))
    }
  )
)