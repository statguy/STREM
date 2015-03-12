library(sp)
library(raster)

StudyArea <- setRefClass(
  Class = "StudyArea",
  fields = list(
    context = "ANY",
    region = "ANY",
    proj4string = "ANY",
    boundary = "SpatialPolygons",
    habitat = "RasterLayer",
    coordinateScale = "integer",
    plotScale = "numeric",
    nBoundarySamples = "integer"
  ),
  methods = list(
    initialize = function(context=NA, region=NA, proj4string=NA, ...) {
      callSuper(context=context, region=region, proj4string=proj4string, ...)
      return(invisible(.self))
    },
    
    setup = function(thin=TRUE, tolerance=0.1, prepareHabitatRaster=FALSE, rawRasterFile) {
      loadBoundary(thin=thin, tolerance=tolerance)
      if (prepareHabitatRaster) {
        if (missing(rawRasterFile)) stop("rawRasterFile parameter missing.")
        prepareHabitatRaster(rawRasterFile=rawRasterFile)
      }
      loadHabitatRaster()
      if (length(plotScale) == 0) plotScale <<- getPlotScale()
      return(invisible(.self))
    },
    
    toGGDF = function() {
      library(ggplot2)
      return(ggplot2::fortify(boundary))
    },
    
    loadBoundary = function(thin, tolerance) {
      stop("Override loadBoundary().")
    },
    
    prepareHabitatRaster = function(rawRasterFile) {
      if (missing(rawRasterFile))
        stop("rawRasterFile or boundary parameter missing.")
      saveHabitatRasterFile(rawRasterFile=rawRasterFile)
      return(invisible(.self))
    },
    
    maskLargeRaster = function(inputRasterFile, outputRasterFile) {
      library(maptools)
      
      message("Writing boundary polygon as shape file...")
      tmpbase <- tempdir()
      tmpfile <- file.path(tmpbase, "boundary")
      writePolyShape(boundary, tmpfile)
      
      ext <- paste(boundary@bbox[1,1], boundary@bbox[2,2], boundary@bbox[1,2], boundary@bbox[2,1])

      message("Cropping raster...")
      tmptif <- tempfile()
      cmd <- paste("gdal_translate -projwin", ext, "-of GTiff", inputRasterFile, tmptif)
      message(cmd)
      status <- system(cmd)
      if (status == 127) stop("GDAL tools not found.")

      message("Masking raster with boundary...")
      tmpshp <- paste(tmpfile, ".shp", sep="")
      cmd <- paste("gdalwarp -co COMPRESS=DEFLATE -co TILED=YES -of GTiff -r lanczos -multi -cutline ", tmpshp, tmptif, outputRasterFile)
      message(cmd)
      status <- system(cmd)
      if (status == 127) stop("GDAL tools not found.")

      #unlink(tmptif)
      #unlink(tmpshp)
      return(invisible(.self))
    },
    
    saveHabitatRasterFile = function(rawRasterFile) {
      if (!file.exists(rawRasterFile))
        stop("Raw habitat raster file ", rawRasterFile, " not found.")
      
      message("Masking raster with boundary...")
      tmpRasterFile <- tempfile()
      maskLargeRaster(rawRasterFile, tmpRasterFile)
      
      message("Storing raster file...")
      file.copy(tmpRasterFile, getHabitatRasterFile(), overwrite=TRUE)
      
      message("Counting habitat frequencies...")
      loadHabitatRaster()
      saveHabitatFrequencies(habitat)
      
      message("Cleaning up...")
      #unlink(tmpRasterFile)
      
      return(invisible(.self))
    },
    
    #getHabitatRasterFile = function() {
    #  stop("Override this method.")
    #},
    
    getHabitatRasterFile = function() {
      return(context$getLongFileName(context$scratchDirectory, name="HabitatRaster", region=region, tag="cropped", ext=".tif"))
    },
    
    loadHabitatRaster = function() {
      library(raster)
      if (!file.exists(getHabitatRasterFile()))
        stop("Habitat raster file ", getHabitatRasterFile(), " not found.")
      habitat <<- raster(getHabitatRasterFile())
      projection(habitat) <<- attr(proj4string, "projargs")
      return(invisible(.self))
    },
            
    findLargestPolygon = function(gadm) {
      library(maptools)
      rownames(gadm@data) <- 1:nrow(gadm@data)      
      p <- list()
      for (i in 1:length(gadm)) {
        largest <- which.max(sapply(gadm@polygons[[i]]@Polygons, slot, "area"))
        p[[i]] <- Polygons(list(Polygon(gadm@polygons[[i]]@Polygons[[largest]])), ID=i)
      }
      p <- SpatialPolygonsDataFrame(SpatialPolygons(p, proj4string=gadm@proj4string), data=gadm@data)
      return(p)
    },
    
    thinGADM = function(gadm, tolerance) {
      library(maptools)
      p <- findLargestPolygon(gadm)
      return(thinnedSpatialPoly(p, tolerance=tolerance))
    },    
    
    loadBoundaryGADM = function(country, level=0, subregion, mainland=FALSE, thin=FALSE, coordinateScale=as.integer(1), tolerance=0.1) {  
      coordinateScale <<- coordinateScale
      
      library(sp)
      library(rgdal)
      
      if (thin==TRUE & mainland==TRUE) mainland <- FALSE
      
      gadm <- try(raster::getData('GADM', country=country, level=level, path=context$scratchDirectory))
      if (inherits(gadm, "try-error"))
        stop("Failed to download country boundary: ", gadm)
      
      if (mainland) gadm <- findLargestPolygon(gadm)
      if (thin) gadm <- thinGADM(gadm, tolerance)
      
      gadm <- spTransform(gadm, proj4string)
      proj4string(gadm) <- proj4string(gadm) # Dunno why is this needed, but it solves some problems
      
      # Scale the coordinates if requested
      if (coordinateScale != 1) {
        sp <- list()
        for (j in seq(along=gadm)) {
          sps <- list()      
          for (i in seq(along=gadm@polygons[[j]]@Polygons))
            sps[[i]] <- Polygon(gadm@polygons[[j]]@Polygons[[i]]@coords * coordinateScale)
          sp[[j]] <- Polygons(sps, ID=gadm@polygons[[j]]@ID)
        }
        sp <- SpatialPolygons(sp, proj4string=gadm@proj4string) # NOTE: Non-standard coordinate scale and proj4 might not work
        gadm <- SpatialPolygonsDataFrame(sp, data=gadm@data)
      }
      
      # Return only the subregion if given
      if (!missing(subregion)) {
        sps <- list()
        for (i in 1:length(subregion))
          sps[[i]] <- Polygons(gadm@polygons[[subregion[i]]]@Polygons, ID=subregion[i])    
        sp <- SpatialPolygons(sps, proj4string=slot(gadm, "proj4string"))
        p <- SpatialPolygonsDataFrame(sp, data=gadm@data[subregion,])
        
        # TODO: remove inner boundaries when joining multiple subregions
        gadm <- p
      }
      
      boundary <<- gadm
      return(invisible(.self))
    },
    
    sampleBoundary = function() {
      if (length(nBoundarySamples) == 0) {
        warning("nBoundarySamples not defined by subclass.")
        return(NULL)
      }
      coords <- coordinates(boundary@polygons[[1]]@Polygons[[1]])
      boundarySL <- SpatialLines(list(Lines(Line(coords), ID=1)), proj4string=proj4string)
      return(spsample(x=boundarySL, n=nBoundarySamples, type="regular"))
    },
    
    getHabitatFrequenciesFileName = function() {
      return(context$getFileName(context$processedDataDirectory, name="HabitatFrequencies", region=studyArea$region))
    },
    
    saveHabitatFrequencies = function(habitatRaster) {
      habitatFrequencies <- freq(habitatRaster, progress="text")
      save(habitatFrequencies, file=getHabitatFrequenciesFileName())
      return(invisible(.self))
    },
    
    loadHabitatFrequencies = function() {
      if (!file.exists(getHabitatFrequenciesFileName()))
        stop("Habitat frequencies file ", getHabitatFrequenciesFileName(), " not found.")
      load(getHabitatFrequenciesFileName())
      return(habitatFrequencies)
    },
    
    getPlotScale = function() {
      return(10 ^ (nchar(as.character(max(dim(habitat)))) - 3))
    },
    
    plotHabitatRaster = function() {
      plot(habitat)
      plot(boundary, add=T)
      return(invisible(.self))
    },
    
    readRasterIntoMemory = function() {
      library(raster)
      if (!inMemory(habitat)) {
        message("Reading habitat raster into memory...")
        habitat <<- readAll(habitat)
      }
      return(invisible(.self))
    },
    
    getMesh = function() {
      stop("Unimplemented method.")      
    },
    
    getQuickMesh = function() {
      stop("Unimplemented method.")
    }
  )
)

TestStudyArea <- setRefClass(
  Class = "TestStudyArea",
  contains = "StudyArea",
  methods = list(
    initialize = function(...) {
      library(sp)
      callSuper(region="test", proj4string=CRS("+init=epsg:2393"), nBoundarySamples=as.integer(0), ...)
      return(.self)
    },
    
    loadBoundary = function(thin, tolerance) {
      x <- -5
      y <- -6
      boundary <<- SpatialPolygons(list(Polygons(list(Polygon(matrix(c(325+x,700+y, 325+x,712+y, 338+x,712+y, 338+x,700+y, 325+x,700+y)*10000, ncol=2, byrow=T))), ID=1)), proj4string=proj4string)
      return(invisible(.self))
    },
    
    loadHabitatRaster = function() {
      callSuper()
      library(raster)
      habitat <<- crop(habitat, boundary)
      habitat <<- readAll(habitat)
      habitat[is.na(habitat[])] <<- 44
      return(invisible(.self))
    },
    
    getHabitatRasterFile = function() {
      return(context$getFileName(context$scratchDirectory, name="HabitatRaster", region="Finland", ext=".tif"))
    },

    getQuickMesh = function() {
      getMesh()
    },
    
    getMesh = function() {
      #list(maxEdge=c(.01e6, .2e6), cutOff=.01e6, coordsScale=1e-6)
      list(coordsScale=1e-6, maxEdge=c(.01e6, .02e6), cutOff=.007e6)
    }
  )
)

FinlandStudyArea <- setRefClass(
  Class = "FinlandStudyArea",
  contains = "StudyArea",
  methods = list(
    initialize = function(...) {
      library(sp)
      callSuper(region="Finland", proj4string=CRS("+init=epsg:2393"), nBoundarySamples=as.integer(40), ...)
      return(invisible(.self))
    },
    
    loadBoundary = function(thin=TRUE, tolerance=0.1) {
      loadBoundaryGADM(country="FIN", thin=thin, tolerance=tolerance)
      return(invisible(.self))
    },
    
    #getQuickMesh = function() { 
    #},
    
    getMesh = function() {
      list(maxEdge=c(.1e6, .2e6), cutOff=.05e6, coordsScale=1e-6)
    }
  )
)

RussiaStudyArea <- setRefClass(
  Class = "RussiaStudyArea",
  contains = "StudyArea",
  methods = list(
    initialize = function(...) {
      library(sp)
      #callSuper(region="Russia", proj4string=CRS("+init=epsg:3410"), ...)
      callSuper(region="Russia", proj4string=CRS("+init=epsg:3857"), ...)
      return(invisible(.self))
    },
    
    loadDistricts = function() {
      library(sp)
      library(maptools)
      library(rgdal)
      
      fileName <- file.path(context$rawDataDirectory, "russia", "adm6_district")
      districts <- readShapeSpatial(fileName)
      proj4string(districts) <- CRS("+proj=longlat +datum=WGS84")
      districts <- spTransform(districts, .self$proj4string)
      
      return(districts)
    },
    
    getBoundaryFileName = function() {
      context$getFileName(dir=context$processedDataDirectory, name="Boundary", region=region)
    },
    
    saveBoundary = function(fileName=getBoundaryFileName(), study) {
      library(sp)
      i <- RussiaWTCIntersections$new(study=study)
      i$loadIntersections()
      districts <- loadDistricts()
      boundary <<- subset(districts, NAME_LAT %in% i$intersections$District_Lat & ADM4_ID %in% i$intersections$RegionID)
      boundary$NAME_ENGLISH <<- boundary$NAME_LAT
      save(boundary, file=fileName)
      return(invisible(.self))
    },
    
    loadBoundary = function(thin=TRUE, tolerance=0.1) {
      library(sp)
      #load(file=getBoundaryFileName(), envir=as.environment(.self))
      loadBoundaryGADM(country="RUS", thin=thin, tolerance=tolerance)
      return(invisible(.self))
    },
    
    loadHabitatRaster = function() {
      library(raster)
      habitat <<- raster(extent(c(2.2,7.5,7.4,11.2)*1e6), nrows=400*100, ncols=1000*100, crs=proj4string)
      return(invisible(.self))
    }
  )
)

FinlandRussiaStudyArea <- setRefClass(
  Class = "FinlandRussiaStudyArea",
  contains = "StudyArea",
  methods = list(
    initialize = function(...) {
      library(sp)
      #callSuper(region="FinlandRussia", proj4string=CRS("+init=epsg:3410"), ...)
      callSuper(region="FinlandRussia", proj4string=CRS("+init=epsg:3857"), ...)
      return(invisible(.self))
    },
    
    loadBoundary = function(thin=TRUE, tolerance=0.1) {
      library(sp)
      library(rgdal)
      
      finland <- FinlandStudyArea$new(context=context)$setup()
      russia <- RussiaStudyArea$new(context=context)$setup()
      
      x <- spTransform(finland$boundary, russia$proj4string)
      y <- russia$boundary
      
      x@data <- data.frame(district=x$NAME_ENGLISH)
      y@data <- data.frame(district=y$NAME_ENGLISH)
      z <- spChFIDs(y, as.character(1:length(y)+1))
      boundary <<- rbind(x, z)

      return(invisible(.self))
    },
    
    loadHabitatRaster = function() {
      library(raster)
      habitat <<- raster(extent(c(2.2,7.5,7.4,11.2)*1e6), nrows=400*100, ncols=1000*100, crs=proj4string)
      return(invisible(.self))
    }
  )
)
