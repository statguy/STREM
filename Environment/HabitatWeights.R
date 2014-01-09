HabitatWeights <- setRefClass(
  Class = "HabitatWeights",
  fields = list(
  ),
  methods = list(
    initialize = function(habitatWeightList) {
    },
    
    classify = function(habitatValues) {
      stop("classify() undefined.")
    },
    
    getWeights = function(habitatValues) {
      stop("getWeights() undefined.")
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
    initialize = function(habitatWeightList=list(urban=1, agriculture=1, forestland=1, peatland=1, water=1)) {
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
          rep(habitatWeightList$urban,13),
          rep(habitatWeightList$agriculture,4),
          rep(habitatWeightList$forestland,18),
          rep(habitatWeightList$peatland,6),
          rep(habitatWeightList$water,3),
          rep(0,256-(1+13+4+18+6+3)) # undefined
        )
      )
    },

    classify = function(habitatValues) {
      x <- weights$type[match(habitatValues, weights$habitat)]
      x[is.na(x)] <- 0
      return(x)
    },

    getWeights = function(habitatValues) {
      x <- weights$weight[match(habitatValues, weights$habitat)]
      x[is.na(x)] <- 0
      return(x)
    }
  )
)
