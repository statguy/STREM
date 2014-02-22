ModelCollection <- setRefClass(
  Class = "ModelCollection",
  fields = list(
    study = "Study",
    name = "character",
    directory = "character",
    modelsList = "list"
  ),
  methods = list(
    initialize = function(...) {
      callSuper(name="Model", ...)
      return(invisible(.self))
    },
    
    getNumberOfModels = function() return(length(modelsList)),
    
    addModel = function(model) {
      modelsList[[getNumberOfModels() + 1]] <<- model
      return(invisible(.self))
    },
    
    getModel = function(index) {
      if (getNumberOfModels() < index)
        stop("Invalid index = ", index, ".")
      return(modelsList[[index]])
    },
    
    applyModels = function(fun, ...) {
      return(lapply(X=modelsList, FUN=fun, ...))
    },
    
    getModelFileName = function(iteration) {
      if (missing(iteration))
        stop("Provide iteration argument.")
      if (inherits(study, "uninitializedField") | length(directory) == 0)
        stop("Set directory and study fields.")

      fileName <- study$context$getLongFileName(dir=directory, name=name, response=study$response, tag=iteration, region=study$studyArea$region)
      return(fileName)
    },
    
    getModelFileNames = function() {
      if (inherits(study, "uninitializedField") | length(directory) == 0)
        stop("Set directory and study fields.")
      
      modelFiles <- study$context$listLongFiles(dir=directory, name=name, response=study$response, tag="\\d+", region=study$studyArea$region)
      if (length(modelFiles) == 0)
        stop("Cannot find any model files in ", directory)
      return(modelFiles)
    },
    
    loadModels = function() {
      stop("Unimplemented method.")
    }
  )
)

SmoothModelCollection <- setRefClass(
  Class = "SmoothModelCollection",
  contains = "ModelCollection",
  fields = list(
  ),
  methods = list(
    initialize = function(...) {
      callSuper(name="SmoothModel", ...)
      return(invisible(.self))
    },
    
    loadModels = function() {
      modelFiles <- getModelFileNames()
      for (iteration in 1:length(modelFiles)) {
        model <- SmoothModel$new(study=study)
        model$loadResult(fileName=modelFiles[iteration])
        addModel(model)
      }
      return(invisible(.self))
    },
    
    collectResults = function(quick=TRUE) {
      applyModels(function(x, quick) x$collectResults(quick=quick), quick=quick)
      return(invisible(.self))
    }
  )
)
