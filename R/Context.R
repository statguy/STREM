Context <- setRefClass(
  "Context",
  fields = list(
    resultDataDirectory = "ANY",
    processedDataDirectory = "ANY",
    rawDataDirectory = "ANY",
    scratchDirectory = "ANY",
    figuresDirectory = "ANY",
    verbose = "logical"
  ),
  methods = list(
    initialize = function(resultDataDirectory=NA, processedDataDirectory=NA, rawDataDirectory=NA, scratchDirectory=NA, figuresDirectory=NA, verbose=FALSE) {
      callSuper(resultDataDirectory=resultDataDirectory, processedDataDirectory=processedDataDirectory, rawDataDirectory=rawDataDirectory, scratchDirectory=scratchDirectory, figuresDirectory=figuresDirectory, verbose=verbose)
    },
    
    getFileName = function(dir, name, response, region, ext=".RData") {
      if (is.na(dir)) stop("Directory has not been specified.")
      
      fileName <- if (missing(response))
        file.path(dir, paste(name, "-", region, ext, sep=""))
      else
        file.path(dir, paste(name, "-", response, "-", region, ext, sep=""))
      if (verbose) message("File = ", fileName)
      return(fileName)
    },
    
    listFiles = function(dir, name, response=NULL, region, tag=NULL, ext=".RData") {
      if (is.na(dir)) stop("Directory has not been specified.")
      pattern <- concat(concat(name, response, region, tag, sep="-"), ext, sep="")
      message("Listing files ", file.path(dir, pattern))
      return(list.files(dir, pattern, full.names=TRUE))
    },
    
    getLongFileName = function(dir, name, response, region, tag, ext=".RData") {
      return(getFileName(dir=dir, name=name, response=response, region=paste(region, tag, sep="-"), ext=ext))
    },
    
    # deprecated
    listLongFiles = function(dir, name, response, region, tag, ext=".RData") {
      return(listFiles(dir=dir, name=name, response=response, region=paste(region, tag, sep="-"), ext=ext))
    },
    
    getIterationIds = function(dir, name, response, region, tag="(\\d+)", ext=".RData") {
      files <- listFiles(dir=dir, name=name, response=response, region=region, tag=tag, ext=ext)
      pattern <- file.path(path.expand(dir), concat(concat(name, response, region, tag, sep="-"), ext, sep=""))
      message("Extraction pattern ", pattern)
      return(as.integer(gsub(pattern, "\\1", files)))
    }
  )
)
