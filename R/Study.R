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
      return(invisible(.self))
    },
    
    newInstance = function(context, isTest=F) {
      context <<- context
      studyArea <<- if (isTest) TestStudyArea$new(context=context)$newInstance()
      else FinlandStudyArea$new(context=context)$newInstance()
      return(invisible(.self))
    },
    
    loadTracksCollection = function() {
      tracks <- SimulatedTracksCollection$new(study=.self)
      tracks$loadTracks()
      return(tracks)
    },
    
    loadIntersectionsCollection = function() {
      intersections <- SimulatedIntersectionsCollection$new(study=.self)
      intersections$loadIntersections()
      return(intersections)      
    },
    
    loadSurveyRoutes = function(n=800) {
      return(FinlandRandomWTCSurveyRoutes$new(study=study)$newInstance(n=800))
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
      return(invisible(.self))
    },
    
    preprocessResponse = function(response, fmiApiKey) {
      library(CNPCluster)
      
      cnpClusterStartLocal()
      
      response <<- response
      
      intersections <- FinlandWTCIntersections$new(study=.self)
      tracks <- FinlandWTCTracks$new(study=.self)
      habitatWeights <- CORINEHabitatWeights$new(study=.self)

      intersections$saveIntersections()
      intersections$saveCovariates(fmiApiKey=fmiApiKey)
      
      tracks$saveTracks()
      
      habitatSelection <- tracks$getHabitatPreferences(habitatWeightsTemplate=habitatWeights, nSamples=30, save=T)
      habitatWeights <- CORINEHabitatWeights$new(study=study)$setHabitatSelectionWeights(habitatSelection)
      habitatWeights$getWeightsRaster(aggregationScale=100, save=T)
      
      cnpClusterStopLocal()
    },
    
    preprocess = function(fmiApiKey) {
      preprocessResponse("canis.lupus", fmiApiKey)
      preprocessResponse("lynx.lynx", fmiApiKey)
      preprocessResponse("rangifer.tarandus.fennicus", fmiApiKey)
      return(invisible(.self))
    },
    
    loadTracks = function() {
      return(FinlandWTCTracks$new(study=.self)$loadTracks())
    },
    
    loadSurveyRoutes = function() {
      return(FinlandWTCSurveyRoutes$new(study=.self)$newInstance())
    },
    
    loadIntersections = function() {
      intersections <- FinlandWTCIntersections$new(study=.self)$loadIntersections()
      intersections$loadCovariates()
      return(intersections)
    },
    
    estimate = function() {
      # TODO
    },
    
    postprocess = function() {
      # TODO
    }
  )
)
