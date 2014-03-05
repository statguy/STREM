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
    
    animate = function(delay=100, name, boundary=study$studyArea$boundary, sameScale=TRUE) {
      library(raster)

      vmin <- min(minValue(rasterStack))
      vmax <- max(maxValue(rasterStack))
      
      #tmpDir <- tempdir()
      for (i in 1:nlayers(rasterStack)) {
        # TODO: same legend scale
        
        layer <- rasterStack[[i]]
        layerName <- names(layer)
        #layerFileName <- file.path(tmpDir, paste("animate", i, ".png", sep=""))
        layerFileName <- context$getFileName(dir=study$context$figuresDirectory, name=paste(name, layerName, sep="-"), response=study$response, region=study$studyArea$region, ext=".png")        
        message("Saving ", layerFileName, "...")
        png(filename=layerFileName, width=dim(layer)[2], height=dim(layer)[1])
        
        if (sameScale) {
          plot(layer, main=layerName,
               col=terrain.colors(99),
               breaks=seq(vmin, vmax, length.out=100))
        }
        else {
          plot(layer, main=layerName, col=terrain.colors(99))
        }
        
        plot(boundary, add=T)
        dev.off()
      }
      
      message("Converting...")
      outputFile <- context$getFileName(dir=study$context$figuresDirectory, name=name, response=study$response, region=study$studyArea$region, ext=".gif")
      #file.path(study$context$figuresDirectory, fileName)
      layerFileNameMask <- context$getFileName(dir=study$context$figuresDirectory, name=paste(name, "*", sep="-"), response=study$response, region=study$studyArea$region, ext=".png")      
      cmd <- paste("convert -loop 0 -delay ", delay, " ", layerFileNameMask, " ", outputFile, sep="")
      #cmd <- paste("convert -loop 0 -delay ", delay, " ", file.path(tmpDir, "animate*.png"), " ", outputFile, sep="")
      message(cmd)
      system(cmd)
      
      return(invisible(.self))
    },
    
    weight = function(weights) {
      for (i in 1:nlayers(rasterStack))
        rasterStack[[i]] <<- rasterStack[[i]] * weights
      return(invisible(.self))
    },
    
    integrate = function(volume=PopulationSize$new(study=study), weights=1) {
      library(raster)
      
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
