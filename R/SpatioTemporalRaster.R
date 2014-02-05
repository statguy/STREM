library(raster)

SpatioTemporalRaster <- setRefClass(
  Class = "SpatioTemporalRaster",
  fields = list(
    meanRaster = "RasterLayer",
    sdRaster = "RasterLayer"
  ),
  methods = list(
    initialize = function(meanRasterList, sdRasterList, ...) {
      callSuper(...)
      if (!missing(meanRasterList) & !missing(sdRasterList))
        setRasters(meanRasterList=meanRasterList, sdRasterList=sdRasterList)
      return(.self)
    },
        
    setRasters = function(meanRasterList, sdRasterList) {
      library(raster)
      meanRaster <<- stack(meanRasterList)
      sdRaster <<- stack(sdRasterList)
    },
    
    addRasters = function(meanRasterLayer, sdRasterLayer) { # TODO: problem here. fix!
      library(raster)
      meanRaster <<- addLayer(meanRaster, meanRasterLayer)
      sdRaster <<- addLayer(sdRaster, sdRasterLayer)
    },
    
    integrate = function(weights=1) {
      library(raster)
      
      volume <- data.frame()
      for (yearIndex in 1:nlayers(meanRaster)) {
        weightedDensity <- meanRaster[[yearIndex]] * weights
        year <- names(meanRaster[[i]])
        volume <- rbind(volume, data.frame(Year=year, Estimated=cellStats(weightedDensity, sum)))
      }
      
      return(invisible(volume))
    }   
  )
)