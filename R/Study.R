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
      tracks$loadTracks()
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
    
    
    preprocessResponse = function(response) {
      response <<- response
      
      intersections <- FinlandWTCIntersections$new(study=.self)
      tracks <- FinlandWTCTracks$new(study=.self)      
      habitatWeights <- CORINEHabitatWeights$new(study=.self)

      intersections$saveIntersections()
      tracks$saveTracks()
      habitatSelection <- tracks$getHabitatPreferences(habitatWeightsTemplate=habitatWeights, nSamples=30, save=T)
      habitatWeights <- CORINEHabitatWeights$new(study=study)$setHabitatSelectionWeights(habitatSelection)
      habitatWeights$getWeightsRaster(aggregationScale=100, save=T)
    },
    
    preprocess = function() {
      preprocessResponse("canis.lupus")
      preprocessResponse("lynx.lynx")
      preprocessResponse("rangifer.tarandus.fennicus")
    }
  )
)
