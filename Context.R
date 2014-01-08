Context <- setRefClass(
  "Context",
  fields = list(
    processedDataDirectory = "ANY",
    scratchDirectory = "ANY"
  ),
  methods = list(
    initialize = function(processedDataDirectory=NA, scratchDirectory=NA) {
      callSuper(processedDataDirectory=processedDataDirectory, scratchDirectory=scratchDirectory)
    }
  )
)
