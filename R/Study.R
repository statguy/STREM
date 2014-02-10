Study <- setRefClass(
  Class = "Study",
  fields = list(
    context = "Context",
    response = "character",
    studyArea = "StudyArea"
  ),
  methods = list(
  )
)

SimulationStudy <- setRefClass(
  Class = "SimulationStudy",
  contains = "Study",
  methods = list(
    initialize = function(...) {
      callSuper(...)
      return(.self)
    },
    
    newInstance = function(context, isTest=F) {
      context <<- context
      studyArea <<- if (isTest) TestStudyArea$new(context=context)$newInstance()
      else FinlandStudyArea$new(context=context)$newInstance()
      return(.self)
    },
    
    loadTracksCollection = function() {
      tracks <- SimulatedTracksCollection$new(study=.self)
      tracks$load()
      return(tracks)
    },
    
    loadIntersectionsCollection = function() {
      intersections <- SimulatedIntersectionsCollection$new(study=.self)
      intersections$load()
      return(intersections)      
    }
  )
)

FinlandWTCStudy <- setRefClass(
  Class = "FinlandWTCStudy",
  contains = "Study",
  fields = list(
  ),
  methods = list(
    initialize = function(context, ...) {
      callSuper(context=context, ...)
      studyArea <<- FinlandStudyArea$new(context=context)$newInstance()
    },
    
    preprocess = function() {
      intersections <- FinlandWTCIntersections$new(study=.self)
      tracks <- FinlandWTCTracks$new(study=.self)
      response <<- "canis.lupus"
      intersections$saveIntersections()
      tracks$saveTracks()
      response <<- "lynx.lynx"
      intersections$saveIntersections()
      tracks$saveTracks()
      response <<- "rangifer.tarandus.fennicus"
      intersections$saveIntersections()
      tracks$saveTracks()
    }
  )
)
