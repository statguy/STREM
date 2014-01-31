Context <- setRefClass(
  "Context",
  fields = list(
    resultDataDirectory = "character",
    processedDataDirectory = "character",
    rawDataDirectory = "character",
    scratchDirectory = "character",
    figuresDirectory = "character"
  ),
  methods = list(
    initialize = function(resultDataDirectory=".", processedDataDirectory=".", rawDataDirectory=".", scratchDirectory=".", figuresDirectory=".") {
      callSuper(resultDataDirectory=resultDataDirectory, processedDataDirectory=processedDataDirectory, rawDataDirectory=rawDataDirectory, scratchDirectory=scratchDirectory, figuresDirectory=figuresDirectory)
    },
    
    getFileName = function(dir, name, response, region, ext=".RData") {
      fileName <- if (missing(response))
        file.path(dir, paste(name, "-", region, ext, sep=""))
      else
        file.path(dir, paste(name, "-", response, "-", region, ext, sep=""))
      message("File = ", fileName)
      return(fileName)
    },
    
    listFiles = function(dir, name, response, region, ext=".RData") {
      pattern <- if (missing(response))
        paste(name, "-", region, ext, sep="")
      else
        paste(name, "-", response, "-", region, ext, sep="")
      return(list.files(dir, pattern))
    }
  )
)
