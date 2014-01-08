Context <- setRefClass(
  "Context",
  fields = list(
    sourceUrl = "character",
    processedDataDirectory = "ANY",
    scratchDirectory = "ANY"
  ),
  methods = list(
    initialize = function(sourceUrl="https://raw.github.com/statguy/R-Winter-Track-Counts/master/", processedDataDirectory=NA, scratchDirectory=NA) {
      callSuper(sourceUrl=sourceUrl, processedDataDirectory=processedDataDirectory, scratchDirectory=scratchDirectory)
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
    }
  )
)
