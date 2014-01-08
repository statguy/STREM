Context <- setRefClass(
  "Context",
  fields = list(
    processedDataDirectory = "character",
    scratchDirectory = "character"
  ),
  methods = list(
    initialize = function(processedDataDirectory, scratchDirectory) {
      processedDataDirectory <<- processedDataDirectory
      scratchDirectory <<- scracthDirectory
    }
  )
)
