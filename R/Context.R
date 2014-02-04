Context <- setRefClass(
  "Context",
  fields = list(
    resultDataDirectory = "ANY",
    processedDataDirectory = "ANY",
    rawDataDirectory = "ANY",
    scratchDirectory = "ANY",
    figuresDirectory = "ANY"
  ),
  methods = list(
    initialize = function(resultDataDirectory=NA, processedDataDirectory=NA, rawDataDirectory=NA, scratchDirectory=NA, figuresDirectory=NA) {
      callSuper(resultDataDirectory=resultDataDirectory, processedDataDirectory=processedDataDirectory, rawDataDirectory=rawDataDirectory, scratchDirectory=scratchDirectory, figuresDirectory=figuresDirectory)
    },
    
    getFileName = function(dir, name, response, region, ext=".RData") {
      if (is.na(dir)) stop("Directory has not been specified.")
      
      fileName <- if (missing(response))
        file.path(dir, paste(name, "-", region, ext, sep=""))
      else
        file.path(dir, paste(name, "-", response, "-", region, ext, sep=""))
      message("File = ", fileName)
      return(fileName)
    },
    
    listFiles = function(dir, name, response, region, ext=".RData") {
      if (is.na(dir)) stop("Directory has not been specified.")
      
      pattern <- if (missing(response))
        paste(name, "-", region, ext, sep="")
      else
        paste(name, "-", response, "-", region, ext, sep="")
      return(list.files(dir, pattern))
    }
  )
)
