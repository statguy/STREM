SpatioTemporalRaster <- setRefClass(
  Class = "SpatioTemporalRaster",
  fields = list(
    study = "Study",
    rasterStack = "RasterStack",
    width = "numeric",
    ext = "character"
  ),
  methods = list(
    initialize = function(layerList, width=16, ext="png", ...) {
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
    
    getRasterFileName = function(name, layerName, ext=.self$ext) {
      return(context$getFileName(dir=study$context$figuresDirectory, name=paste(name, layerName, sep="-"), response=study$response, region=study$studyArea$region, ext=paste(".", ext, sep="")))
    },
    
    saveRasterFile = function(p, layer, name, layerName, ext=.self$ext, ...) {
      if (missing(p) | missing(name) | missing(layerName))
        stop("Argument p, name or layerName missing.")
      aspectRatio <- dim(layer)[1] / dim(layer)[2]
      height <- 12
      if (aspectRatio > 1) width0 <- height / aspectRatio
      else height <- width * aspectRatio
      ggsave(p, filename=getRasterFileName(name=name, layerName=layerName, ext=ext), width=width0, height=height, ...)
    },
    
    plotLayer = function(layerName, plotTitle, legendTitle, digits=2, breaks=round(seq(minValue(rasterStack[[layerName]]), maxValue(rasterStack[[layerName]]), length.out=7), digits=digits), boundary=FALSE, plot=TRUE, save=FALSE, name, ...) {
      library(ggplot2)
      library(raster)
      library(rasterVis)

      message("Breaks:")
      print(breaks)
      
      layer <- crop(rasterStack[[layerName]], extent(study$studyArea$boundary))
      p <- gplot(layer) + geom_raster(aes(fill=value)) + coord_equal() +
        #scale_fill_gradientn(colours=terrain.colors(99), na.value=NA)
        scale_fill_gradient(low="white", high="steelblue", na.value=NA, breaks=breaks)
      p <- if (!missing(legendTitle)) p + theme_raster(legend.position="right") + guides(fill=guide_legend(title=legendTitle))
      else p + theme_raster()
      
      if (!missing(boundary)) {
        boundaryDF <- if (class(boundary) == "logical") ggplot2::fortify(study$studyArea$boundary)
        else boundary
        p <- p + geom_polygon(data=boundaryDF, aes(long, lat), colour="black", fill=NA)
      }
      
      if (!missing(plotTitle)) p <- p + ggtitle(plotTitle)
      
      if (plot) print(p)
      if (save) saveRasterFile(p=p, layer=layer, name=name, layerName=layerName, ...)
      
      return(invisible(p))
    },
    
    interpolate = function(xyzt, timeVariable=colnames(xyzt)[4], templateRaster=study$getTemplateRaster(), transform=identity, inverseTransform=identity, layerNames, boundary, weights=1) {
      library(ST)
      library(raster)
      library(plyr)
      
      rasterStack <<- ST::multiRasterInterpolate(xyzt, variables=timeVariable, templateRaster=templateRaster, transform=transform, inverseTransform=inverseTransform)
      if (!missing(boundary)) {
        rasterStack <<- stack(llply(1:nlayers(rasterStack),
                              function(i, rasterStack, boundary) {
                                message("Cropping ", i, "/", nlayers(rasterStack), "...")
                                mask(rasterStack[[i]], boundary)
                              },
                              rasterStack=rasterStack, boundary=boundary, .parallel=T))
        #for (i in 1:nlayers(rasterStack)) {          
        #  rasterStack[[i]] <<- mask(rasterStack[[i]], boundary)
        #}
      }
      weight(weights)
      if (!missing(layerNames)) names(rasterStack) <<- layerNames
      return(invisible(.self))
    },
    
    getColorBreaks = function(n=6) {
      vmin <- min(minValue(rasterStack))
      vmax <- max(maxValue(rasterStack))
      return(seq(vmin, vmax, length.out=n))
    },
    
    getColorScale = function(n=6) {
      return(terrain.colors(n-1))
    },
    
    animate = function(name, delay=100, boundary=FALSE, sameScale=TRUE, nColor=6, legend, convertCmd="convert", ggfun, ggfunParams, ...) {
      if (missing(name)) stop("Argument name missing.")
      
      library(ggplot2)
      library(raster)
      library(rasterVis)
      
      if (boundary) boundaryDF <- fortify(study$studyArea$boundary)
      colors <- getColorScale(n=nColor)
      breaks <- getColorBreaks(n=nColor)
      
      for (i in 1:nlayers(rasterStack)) {
        layer <- crop(rasterStack[[i]], extent(study$studyArea$boundary))
        layerName <- names(layer)
        if (substr(layerName, start=1, stop=1) == "X") layerName <- substr(layerName, start=2, stop=nchar(layerName))

        message("Processing ", layerName, ": range = ", minValue(layer), "/", min(minValue(rasterStack)), " - ", maxValue(layer), "/", max(maxValue(rasterStack)))
        p <- gplot(layer) + geom_raster(aes(fill=value)) + coord_equal() + ggtitle(layerName)
        p <- if (!missing(legend)) { legendTitle <- legend; p + theme_raster(legend.position="right") }
        else { legendTitle="Density"; p + theme_raster() } 
        p <- if (sameScale) p + scale_fill_gradientn(legendTitle, colours=colors, limits=range(breaks), na.value=NA)
        else p + scale_fill_gradientn(legendTitle, colours=colors, na.value=NA)
        if (boundary) p <- p + geom_polygon(data=boundaryDF, aes(long, lat), colour="white", fill=NA)
        if (!missing(ggfun)) {
          ggfun <- match.fun(ggfun)
          p <- p + ggfun(i=i, params=ggfunParams)
        }
        print(p)

        saveRasterFile(p=p, layer=layer, name=name, layerName=layerName, ...)
        if (ext != "png") saveRasterFile(p=p, layer=layer, name=name, layerName=layerName, ext="png", ...)
      }
      
      message("Converting...")
      outputFile <- context$getFileName(dir=study$context$figuresDirectory, name=name, response=study$response, region=study$studyArea$region, ext=".gif")
      layerFileNameMask <- getRasterFileName(name=name, layerName="*", ext="png")
      cmd <- paste(convertCmd, " -loop 0 -delay ", delay, " ", layerFileNameMask, " ", outputFile, sep="")
      message(cmd)
      x <- system(cmd)
      if (x > 0) stop("System command failed.")
      
      return(invisible(.self))
    },
    
    weight = function(weights) {
      for (i in 1:nlayers(rasterStack)) {
        name <- names(rasterStack[[i]])
        if (inherits(weights, "RasterStack") & length(weights) > 1) rasterStack[[i]] <<- rasterStack[[i]] * weights[[i]]
        else rasterStack[[i]] <<- rasterStack[[i]] * weights
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
