SurveyRoutes <- setRefClass(
  Class = "SurveyRoutes",
  fields = list(
    study = "Study",
    surveyRoutes = "SpatialLines",
    centroids = "SpatialPoints",
    lengths = "numeric",
    lengthsByHabitat = "list",
    DEBUGinnercellLengths = "list"
  ),
  methods = list(    
    plotSurveyRoutes = function() {
      plot(study$studyArea$boundary)
      plot(surveyRoutes, col="blue", add=T)
      return(invisible(.self))
    },
    
    toGGDF = function() {
      library(sp)
      library(ggplot2)
      surveyRoutesSPDF <- SpatialLinesDataFrame(surveyRoutes, data=data.frame(dummy=1:length(surveyRoutes)), match.ID=FALSE)
      return(ggplot2::fortify(surveyRoutesSPDF))
    },
    
    getSurveyRoutesFileName = function() {
      return(study$context$getFileName(dir=study$context$resultDataDirectory, name="SurveyRoutes", response=study$response, region=study$studyArea$region))
    },
    
    loadSurveyRoutes = function(fileName=getSurveyRoutesFileName(), findLengths=FALSE) {
      load(fileName, env=as.environment(.self))
      if (findLengths) getLengths()
      return(invisible(.self))
    },

    saveSurveyRoutes = function(fileName=getSurveyRoutesFileName(), findLengths=TRUE) {
      if (findLengths) getLengths()
      save(surveyRoutes, centroids, lengths, lengthsByHabitat, file=fileName)
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
    
    cutSurveyRoutesByBoundary = function() {
      cutSurveyRoutes <- gIntersection(study$studyArea$boundary, surveyRoutes, byid=T, id=names(surveyRoutes))
      return(cutSurveyRoutes)
    },
    
    getLengths = function() {
      cutSurveyRoutes <- cutSurveyRoutesByBoundary()
      lengths <<- SpatialLinesLengths(cutSurveyRoutes)
      return(invisible(.self))
    },
    
    getlengthsByHabitatType = function(habitatWeights, readHabitatIntoMemory=TRUE, debug=FALSE) {
      if (!inherits(habitatWeights, "HabitatWeights"))
        stop("Parameter habitatClassification must be class of HabitatWeights")
      
      if (readHabitatIntoMemory) study$studyArea$readRasterIntoMemory()
      
      cutSurveyRoutes <- cutSurveyRoutesByBoundary()
      innercellLengths <- findInnercellLengths(study$studyArea$habitat, cutSurveyRoutes)
      if (debug==TRUE) DEBUGinnercellLengths <<- innercellLengths
      
      lengthsByHabitat <<- lapply(innercellLengths, function(x, habitatWeights) {
        classes <- habitatWeights$classify(x[,2])
        cd <- cbind(classes, x[,3])
        hl <- aggregate(cd[,2], sum, by=list(cd[,1]))
        missingIndex <- !habitatWeights$habitatTypes %in% hl[,1]
        hl.final <- data.frame()
        j <- 1
        for (i in 1:length(missingIndex)) {
          if (missingIndex[i] == FALSE) {
            hl.final <- rbind(hl.final, hl[j,])
            j <- j + 1
          }
          else hl.final <- rbind(hl.final, data.frame(Group.1=i, x=0))
        }
        
        hl.final[,2]
      }, habitatWeights=habitatWeights)
      
      fullLengths <- sapply(lengthsByHabitat, function(x) sum(x))
      if (any(fullLengths != lengths))
        warning("Lengths mismatch.")
      
      #weightedEffort <- sapply(lengthsByHabitat, function(x) sum(x * habitatWeights$getHabitatSelectionWeights()))
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
      candidateArea <- RandomInitialPopulation$new(studyArea=study$studyArea)
      centroids <<- candidateArea$randomize(nSurveyRoutes)
      angles <- runif(length(centroids), 0, 2*pi)
      surveyRoutes <<- getTriangles(centroids, angles, 4000)
      getLengths()
      if (save) saveSurveyRoutes()
      return(invisible(.self))
    },
    
    getSurveyRoutesFileName = function() {
      return(study$context$getFileName(dir=study$context$resultDataDirectory, name="SurveyRoutes", response="FinlandRandomWTC", region=study$studyArea$region))
    }
  )
)

FinlandRandomForestWTCSurveyRoutes <- setRefClass(
  Class = "FinlandRandomForestWTCSurveyRoutes",
  contains = "SurveyRoutes",
  fields = list(
  ),
  methods = list(
    initialize = function(...) {
      callSuper(...)
      return(invisible(.self))
    },
    
    randomizeSurveyRoutes = function(nSurveyRoutes, save=FALSE, isTest=FALSE) {
      library(plyr)
      library(sp)
      library(raster)
      library(rgeos)
      
      candidateArea <- RandomInitialPopulation$new(studyArea=study$studyArea)
      
      triangles <- llply(1:nSurveyRoutes, function(i, nSurveyRoutes, candidateArea, habitat) {
        message("Randomizing survey route ", i, "/", nSurveyRoutes, "...")
        
        habitatWeights <- CORINEHabitatWeights$new()
        candidateCentroids <- candidateArea$randomize(100) # Max tries
        angles <- runif(length(candidateCentroids), 0, 2*pi)
        candidateSurveyRoutes <- getTriangles(candidateCentroids, angles, 4000)
        for (j in 1:length(candidateSurveyRoutes)) {
          habitatTypes <- raster::extract(habitat, candidateSurveyRoutes[j])[[1]]
          classifiedHabitatTypes <- habitatWeights$classify(habitatTypes) == 3 # Forest
          if (sum(classifiedHabitatTypes) / length(classifiedHabitatTypes) > .9) {
            triangle <- candidateSurveyRoutes[j]@lines[[1]]
            triangle@ID <- as.character(i)
            return(triangle)
          }
        }
        
        stop("Failed to find survey routes with the given condition.")
      }, nSurveyRoutes=nSurveyRoutes, candidateArea=candidateArea, habitat=study$studyArea$habitat, .parallel=T)
      
      surveyRoutes <<- SpatialLines(triangles, proj4string=study$studyArea$proj4string)
      centroids <<- gCentroid(surveyRoutes, byid=TRUE)
      getLengths()
      if (save) saveSurveyRoutes()
      return(invisible(.self))
    },
    
    getSurveyRoutesFileName = function() {
      return(study$context$getFileName(dir=study$context$resultDataDirectory, name="SurveyRoutes", response="FinlandRandomForestWTC", region=study$studyArea$region))
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
    
    loadSurveyRoutes = function(context, fileName=getSurveyRoutesFileName(), findLengths=TRUE, nSurveyRoutes) {
      study0 <- FinlandWTCStudy$new(context=context, response="canis.lupus")
      intersections <- FinlandWTCIntersections$new(study=study0)$loadIntersections()
      centroids <<- intersections$getSurveyLocations()
      # Angles are unknown, randomize. TODO: Find the most likely angles given the landscape.
      angles <- runif(length(centroids), 0, 2 * pi)
      surveyRoutes <<- getTriangles(centroids, angles, 4000)
      if (findLengths) {
        getLengths()
        index <- lengths > 0
        surveyRoutes <<- surveyRoutes[index]
        centroids <<- centroids[index]
        lengths <<- lengths[index]
      }
        
      if (!missing(nSurveyRoutes)) {
        index <- sample(1:length(surveyRoutes), nSurveyRoutes)
        #index <- 1:nSurveyRoutes
        surveyRoutes <<- surveyRoutes[index]
        centroids <<- centroids[index]
        lengths <<- lengths[index]
      }
      
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
      return(study$context$getFileName(dir=study$context$resultDataDirectory, name="SurveyRoutes", response="Random", region=study$studyArea$region))
    }
  )
)