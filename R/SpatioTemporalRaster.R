library(raster)

SpatioTemporalRaster <- setRefClass(
  Class = "SpatioTemporalRaster",
  fields = list(
    rasterStack = "RasterStack"
  ),
  methods = list(
    initialize = function(meanRasterList, ...) {
      callSuper(...)
      if (!missing(meanRasterList))
        setLayers(meanRasterList=meanRasterList)
      return(.self)
    },
        
    setLayers = function(layerList) {
      library(raster)
      rasterStack <<- stack(layerList)
    },
    
    addLayer = function(layer) {
      library(raster)
      rasterStack <<- raster::addLayer(rasterStack, layer)
    },
    
    integrate = function(weights=1) {
      library(raster)
      
      volume <- PopulationSize$new()
      
      for (index in 1:nlayers(rasterStack)) {
        weightedDensity <- rasterStack[[index]] * weights
        weightedSum <- cellStats(weightedDensity, sum)
        year <- names(rasterStack[[index]])
        if (substr(year, 1, 1) == "X") year <- substr(year, 2, nchar(year))
        volume$addYearSize(year=year, size=weightedSum)
      }
      
      return(volume)
    }
  )
)