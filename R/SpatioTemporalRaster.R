SpatioTemporalRaster <- setRefClass(
  Class = "SpatioTemporalRaster",
  fields = list(
    study = "Study",
    rasterStack = "RasterStack",
    width = "numeric",
    ext = "character"
  ),
  methods = list(
    initialize = function(layerList, width=8, ext="png", ...) {
      callSuper(...)
      width <<- width
      ext <<- ext
      if (!missing(layerList))
        setLayers(layerList=layerList)
      return(invisible(.self))
    },
        
    setLayers = function(layerList) {
      library(raster)
      rasterStack <<- stack(layerList)
      return(invisible(.self))
    },
    
    addLayer = function(layer, name) {
      library(raster)
      if (!missing(name)) names(layer) <- name
      rasterStack <<- raster::addLayer(rasterStack, layer)
      return(invisible(.self))
    },
    
    getRasterFileName = function(name, layerName) {
      return(context$getFileName(dir=study$context$figuresDirectory, name=paste(name, layerName, sep="-"), response=study$response, region=study$studyArea$region, ext=paste(".", ext, sep="")))
    },
    
    saveRasterFile = function(p, layer, name, layerName) {
      if (missing(p) | missing(name) | missing(layerName))
        stop("Argument p, name or layerName missing.")
      aspectRatio <- dim(layer)[2] / dim(layer)[1]
      ggsave(p, filename=getRasterFileName(name=name, layerName=layerName), width=width, height=width * aspectRatio)
    },
    
    plotLayer = function(layerName, plotTitle, legendTitle, boundary=FALSE, save=FALSE, name) {
      library(ggplot2)
      library(raster)
      library(rasterVis)

      layer <- crop(rasterStack[[layerName]], extent(study$studyArea$boundary))
      p <- gplot(layer) + geom_raster(aes(fill=value)) + coord_equal()
        #scale_fill_gradientn(colours=terrain.colors(99), na.value=NA)
      p <- if (!missing(legendTitle)) p + theme_raster(legend.position="right") + guides(fill=guide_legend(title=legendTitle))
      else p + theme_raster()
      if (boundary) {
        boundaryDF <- ggplot2::fortify(study$studyArea$boundary)
        p <- p + geom_polygon(data=boundaryDF, aes(long, lat), colour="white", fill=NA)
      }
      if (!missing(plotTitle)) p <- p + ggtitle(plotTitle)
      
      print(p)
      if (save) saveRasterFile(p=p, layer=layer, name=name, layerName=layerName)
      
      return(invisible(p))
    },
    
    interpolate = function(xyzt, timeVariable, templateRaster=study$getTemplateRaster(), transform=identity, inverseTransform=identity, layerNames) {
      library(ST)
      if (missing(timeVariable)) {
        if (missing(layerNames)) stop("Missing layerNames argument.")
        x <- rasterInterpolate(xyzt, templateRaster=templateRaster, transform=transform, inverseTransform=inverseTransform)
        addLayer(x, layerNames)
      }
      else {
        rasterStack <<- multiRasterInterpolate(xyzt, variables=timeVariable, templateRaster=templateRaster, transform=transform, inverseTransform=inverseTransform)
        if (!missing(layerNames)) names(rasterStack) <<- layerNames
      }
      return(invisible(.self))
    },
    
    animate = function(name, delay=100, boundary=FALSE, sameScale=TRUE) {
      if (missing(name)) stop("Argument name missing.")
      
      library(raster)
      library(rasterVis)
      
      vmin <- min(minValue(rasterStack))
      vmax <- max(maxValue(rasterStack))
      if (boundary) boundaryDF <- fortify(study$studyArea$boundary)
      
      for (i in 1:nlayers(rasterStack)) {
        layer <- crop(rasterStack[[i]], extent(study$studyArea$boundary))
        layerName <- names(layer)
        if (substr(layerName, start=1, stop=1) == "X") layerName <- substr(layerName, start=2, stop=nchar(layerName))
        #layerFileName <- context$getFileName(dir=study$context$figuresDirectory, name=paste(name, layerName, sep="-"), response=study$response, region=study$studyArea$region, ext=".png")        
        #message("Saving ", layerFileName, "...")
        
        p <- gplot(layer) + geom_raster(aes(fill=value)) + coord_equal() + theme_raster() + ggtitle(layerName)
        p <- if (sameScale) p + scale_fill_gradientn(colours=terrain.colors(99), breaks=seq(vmin, vmax, length.out=100), na.value=NA)
        else p + scale_fill_gradientn(colours=terrain.colors(99), na.value=NA)
        if (boundary) p <- p + geom_polygon(data=boundaryDF, aes(long, lat), colour="white", fill=NA)
        
        print(p)
        #ggsave(p, filename=layerFileName, width=width, height=width * aspectRatio)
        saveRasterFile(p=p, layer=layer, name=name, layerName=layerName)
      }
      
      message("Converting...")
      outputFile <- context$getFileName(dir=study$context$figuresDirectory, name=name, response=study$response, region=study$studyArea$region, ext=".gif")
      layerFileNameMask <- getRasterFileName(name=name, layerName="*")
      #layerFileNameMask <- context$getFileName(dir=study$context$figuresDirectory, name=paste(name, "*", sep="-"), response=study$response, region=study$studyArea$region, ext=".png")      
      cmd <- paste("convert -loop 0 -delay ", delay, " ", layerFileNameMask, " ", outputFile, sep="")
      message(cmd)
      system(cmd)
      
      return(invisible(.self))
    },
    
    weight = function(weights) {
      for (i in 1:nlayers(rasterStack)) {
        name <- names(rasterStack[[i]])
        rasterStack[[i]] <<- rasterStack[[i]] * weights
        names(rasterStack[[i]]) <<- name
      }
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
    },
    
    show = function() {
      print(rasterStack)
      return(invisible(.self))
    }
  )
)
