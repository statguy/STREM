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
    
    addLayer = function(layer) { # TODO: problem here. fix!
      library(raster)
      rasterStack <<- raster::addLayer(rasterStack, layer)
    },
    
    integrate = function(weights=1) {
      library(raster)
      
      volume <- data.frame()
      for (index in 1:nlayers(rasterStack)) {
        weightedDensity <- rasterStack[[index]] * weights
        year <- names(meanRaster[[index]])
        volume <- rbind(volume, data.frame(Year=year, Estimated=cellStats(weightedDensity, sum)))
      }
      
      return(invisible(volume))
    }   
  )
)