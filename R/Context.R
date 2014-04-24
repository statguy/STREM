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
    
    listFiles = function(dir, name, response, region, ext=".RData") {
      if (is.na(dir)) stop("Directory has not been specified.")
      
      pattern <- if (missing(response))
        paste(name, "-", region, ext, sep="")
      else
        paste(name, "-", response, "-", region, ext, sep="")
      return(list.files(dir, pattern, full.names=TRUE))
    },
    
    getLongFileName = function(dir, name, response, region, tag, ext=".RData") {
      return(getFileName(dir=dir, name=name, response=response, region=paste(region, tag, sep="-"), ext=ext))
    },
    
    listLongFiles = function(dir, name, response, region, tag, ext=".RData") {
      return(listFiles(dir=dir, name=name, response=response, region=paste(region, tag, sep="-"), ext=ext))
    },
    
    getIterationIds = function(dir, name, response, region, tag="\\d+", ext=".RData") {
      files <- listLongFiles(dir=dir, name=name, response=response, region=region, tag=tag, ext=ext)
      pattern <- file.path(path.expand(dir), paste(paste(name, response, region, "(\\d+)", sep="-"), ext, sep=""))
      return(as.integer(gsub(pattern, "\\1", files)))
    }
  )
)
