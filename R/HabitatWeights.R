HabitatWeights <- setRefClass(
  Class = "HabitatWeights",
  fields = list(
  ),
  methods = list(
    initialize = function(...) {
      callSuper(...)
      return(.self)
    },
    
    classify = function(habitatValues) {
      return(rep(1, times=length(habitatValues)))
    },
    
    getWeights = function(habitatValues) {
      return(rep(1, times=length(habitatValues)))
    },
    
    getRaster = function(habitat, aggregationScale=100) { # TODO: determine aggregation scale automatically
      library(raster)
      weightsRaster <- aggregate(habitat, aggregationScale,
        function(habitatValue, na.rm) mean(getWeights(habitatValue), na.rm=na.rm), na.rm=T)
      return(weightsRaster)
    }
  )
)

CORINEHabitatWeights <- setRefClass(
  Class = "CORINEHabitatWeights",
  contains = "HabitatWeights",
  fields = list(
    weights = "data.frame"
  ),
  methods = list(
    initialize = function(habitatWeightsList=list(urban=1, agriculture=1, forestland=1, peatland=1, water=1), ...) {
      callSuper(...)
      weights <<- data.frame(
        habitat=0:255,
        type=c(
          rep(0,1),  # unknown
          rep(1,13), # urban
          rep(2,4),  # agriculture
          rep(3,18), # forestland
          rep(4,6),  # peatland
          rep(5,3),  # water
          rep(0,256-(1+13+4+18+6+3)) # undefined
        ),
        weight=c(
          rep(0,1), # unknown
          rep(habitatWeightsList$urban,13),
          rep(habitatWeightsList$agriculture,4),
          rep(habitatWeightsList$forestland,18),
          rep(habitatWeightsList$peatland,6),
          rep(habitatWeightsList$water,3),
          rep(0,256-(1+13+4+18+6+3)) # undefined
        )
      )
    },

    classify = function(habitatValues, na.value=NA) {
      x <- weights$type[match(habitatValues, weights$habitat)]
      x[is.na(x)] <- na.value
      return(x)
    },

    getWeights = function(habitatValues, na.value=NA) {
      x <- weights$weight[match(habitatValues, weights$habitat)]
      x[is.na(x)] <- na.value
      return(x)
    }
  )
)
