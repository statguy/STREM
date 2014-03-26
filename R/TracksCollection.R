TracksCollection <- setRefClass(
  Class = "TracksCollection",
  fields = list(
    study = "Study",
    tracksList = "list"
  ),
  methods = list(
    initialize = function(...) {
      callSuper(...)
      return(invisible(.self))
    },

    addTracks = function(tracks) {
      tracksList[[length(tracksList) + 1]] <<- tracks
      return(invisible(.self))
    },
    
    getTracks = function(index) {
      if (length(tracksList) < index) stop("Invalid index = ", index, ".")
      return(tracksList[[index]])
    },
    
    getTracksFileName = function() {
      return(context$getFileName(context$resultDataDirectory, name="TracksCollection", response=study$response, region=study$studyArea$region))
    },
    
    loadTracks = function(fileName=getTracksFileName()) {
      load(file=fileName, envir=as.environment(.self))
    },
    
    saveTracks = function(fileName=getTracksFileName()) {
      save(tracksList, file=fileName)
    }
  )
)

SimulatedTracksCollection <- setRefClass(
  Class = "SimulatedTracksCollection",
  contains = "TracksCollection",
  methods = list(
  )
)
