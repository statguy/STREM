Context <- setRefClass(
  "Context",
  fields = list(
    sourceUrl = "character",
    resultDataDirectory = "ANY",
    processedDataDirectory = "ANY",
    rawDataDirectory = "ANY",
    scratchDirectory = "ANY",
    figuresDirectory = "ANY"
  ),
  methods = list(
    initialize = function(sourceUrl="https://raw.github.com/statguy/R-Winter-Track-Counts/master/", resultDataDirectory=NA, processedDataDirectory=NA, rawDataDirectory=NA, scratchDirectory=NA, figuresDirectory=NA) {
      callSuper(sourceUrl=sourceUrl, resultDataDirectory=resultDataDirectory, processedDataDirectory=processedDataDirectory, rawDataDirectory=rawDataDirectory, scratchDirectory=scratchDirectory, figuresDirectory=figuresDirectory)
    },
    
    loadSource = function(sourceFile) {
      if (substr(sourceUrl, 1, 4) == "http") {
        require(devtools)
        src <- paste(sourceUrl, sourceFile, sep="/")
        message("Loading source from ", src)
        source_url(src)
      }
      else {
        src <- file.path(sourceUrl, sourceFile)
        message("Loading source from ", src)
        source(src)
      }
    },
    
    getFileName = function(dir, name, response, region, ext=".RData") {
      fileName <- if (missing(response))
        file.path(dir, paste(name, "-", region, ext, sep=""))
      else
        file.path(dir, paste(name, "-", response, "-", region, ext, sep=""))
      message("File = ", fileName)
      return(fileName)
    }
  )
)
