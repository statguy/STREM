Data <- setRefClass(
  Class = "Data",
  fields = list(
    context = "Context",
    study = "Study",
    dataName = "character",
    data = "ANY"
  ),
  methods = list(
    initialize = function(preprocessData=FALSE, ...) {
      callSuper(...)
      if (preprocessData) saveData()
      return(.self)
    },
    
    newInstance = function() {
      data <<- loadData()
      return(.self)
    },
    
    getDataFileName = function() {
      return(context$getFileName(context$processedDataDirectory, name=dataName, response=study$response, region=study$studyArea$region))
    },
    
    saveData = function() {
      stop("Override saveData() method.")
    },
    
    loadData = function() {
      load(getDataFileName())
      return(data)
    }
  )
)
