SpatioTemporalRaster <- setRefClass(
  Class = "SpatioTemporalRaster",
  fields = list(
    study = "Study",
    rasterStack = "RasterStack"
  ),
  methods = list(
    initialize = function(meanRasterList, ...) {
      callSuper(...)
      if (!missing(meanRasterList))
        setLayers(meanRasterList=meanRasterList)
      return(invisible(.self))
    },
        
    setLayers = function(layerList) {
      library(raster)
      rasterStack <<- stack(layerList)
      return(invisible(.self))
    },
    
    addLayer = function(layer) {
      library(raster)
      rasterStack <<- raster::addLayer(rasterStack, layer)
      return(invisible(.self))
    },
    
    interpolate = function(xyzt, timeVariable, templateRaster=study$getTemplateRaster(), transform=identity, inverseTransform=identity) {
      library(ST)
      rasterStack <<- multiRasterInterpolate(xyzt, variables=timeVariable, templateRaster=templateRaster, transform=transform, inverseTransform=inverseTransform)        
      return(invisible(.self))
    },
    
    animate = function(delay=100, fileName, boundary=study$studyArea$boundary) {
      library(raster)

      tmpDir <- tempdir()
      for (i in 1:nlayers(rasterStack)) {
        # TODO: same legend scale
        
        layerFileName <- file.path(tmpDir, paste("animate", i, ".png", sep=""))
        message("Saving ", layerFileName, "...")
        layer <- rasterStack[[i]]
        png(filename=layerFileName, width=dim(layer)[2], height=dim(layer)[1])
        plot(layer, main=names(layer))
        plot(boundary, add=T)
        dev.off()
      }
      
      message("Converting...")
      outputFile <- file.path(study$context$figuresDirectory, fileName)
      cmd <- paste("convert -loop 0 -delay ", delay, " ", file.path(tmpDir, "animate*.png"), " ", outputFile, sep="")
      message(cmd)
      system(cmd)
      
      return(invisible(.self))
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
