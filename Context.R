Context <- setRefClass(
  "Context",
  fields = list(
    processedDataDirectory = "character",
    scratchDirectory = "character"
  ),
  methods = list(
    initialize = function(processedDataDirectory, scratchDirectory) {
      .self$processedDataDirectory = processedDataDirectory
      .self$scratchDirectory = scracthDirectory
    }
  )
)
