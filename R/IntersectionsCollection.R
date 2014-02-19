IntersectionsCollection <- setRefClass(
  Class = "IntersectionsCollection",
  fields = list(
    study = "Study",
    intersectionsList = "list"
  ),
  methods = list(
    initialize = function(...) {
      callSuper(...)
    },
    
    getNumberOfIntersections = function() return(length(intersectionsList)),
    
    addIntersections = function(intersections) {
      intersectionsList[[getNumberOfIntersections() + 1]] <<- intersections
    },
    
    getIntersections = function(index) {
      if (getNumberOfIntersections() < index) stop("Invalid index = ", index, ".")
      return(intersectionsList[[index]])
    },
    
    loadIntersections = function() {
      stop("Unimplemented method.")
    },
    
    applyIntersections = function(fun, ...) {
      return(lapply(X=intersectionsList, FUN=fun, ...))
    },
    
    estimate = function() {
      stop("Unimplemented method.")
    }
  )
)

SimulatedIntersectionsCollection <- setRefClass(
  Class = "SimulatedIntersectionsCollection",
  contains = "IntersectionsCollection",
  methods = list(
    loadIntersections = function() {     
      intersectionFiles <- context$listLongFiles(dir=study$context$resultDataDirectory, name="Intersections", response=study$response, tag="\\d+", region=study$studyArea$region)
      
      for (iteration in 1:length(intersectionFiles)) { # TODO: better to use iteration number rather than number of files
        intersections <- SimulatedIntersections$new(study=study, iteration=iteration)
        intersections$loadIntersections()
        addIntersections(intersections)
      }
    },
    
    estimate = function(meshParams, save=FALSE) {
      models <- ModelCollection$new(study=study, directory=study$context$scratchDirectory)
      
      for (iteration in 1:getNumberOfIntersections()) {
        model <- SmoothModel$new(study=study)$setup(intersections=intersectionsList[[iteration]], meshParams=meshParams)
        #model$plotMesh(surveyRoutes=surveyRoutes)
        model$estimate(save=save, fileName=models$getModelFileName(iteration))
        models$addModel(model)
      }
      
      return(invisible(models))
    }
  )
)
