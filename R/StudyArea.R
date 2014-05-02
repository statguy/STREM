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
    plotScale = "numeric"
  ),
  methods = list(
    initialize = function(context=NA, region=NA, proj4string=NA, ...) {
      callSuper(context=context, region=region, proj4string=proj4string, ...)
      return(invisible(.self))
    },
    
    newInstance = function(thin=TRUE, tolerance=0.1, prepareHabitatRaster=FALSE, rawRasterFile) {
      loadBoundary(thin=thin, tolerance=tolerance)
      if (prepareHabitatRaster) {
        if (missing(rawRasterFile)) stop("rawRasterFile parameter missing.")
        prepareHabitatRaster(rawRasterFile=rawRasterFile)
      }
      loadHabitatRaster()
      if (length(plotScale) == 0) plotScale <<- getPlotScale()
      return(invisible(.self))
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
      return(context$getFileName(context$scratchDirectory, name="HabitatRaster", region=region, ext=".tif"))
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
      if (!inMemory(habitat)) {
        message("Reading habitat raster into memory...")
        habitat <<- readAll(habitat)
      }
      return(invisible(.self))
    }    
  )
)

TestStudyArea <- setRefClass(
  Class = "TestStudyArea",
  contains = "StudyArea",
  methods = list(
    initialize = function(...) {
      callSuper(region="test", proj4string=CRS("+init=epsg:2393"), ...)
      return(.self)
    },
    
    loadBoundary = function(thin, tolerance) {
      x <- -5
      boundary <<- SpatialPolygons(list(Polygons(list(Polygon(matrix(c(325+x,700, 325+x,712, 338+x,712, 338+x,700, 325+x,700)*10000, ncol=2, byrow=T))), ID=1)), proj4string=proj4string)
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
    }
  )
)

FinlandStudyArea <- setRefClass(
  Class = "FinlandStudyArea",
  contains = "StudyArea",
  methods = list(
    initialize = function(...) {
      callSuper(region="Finland", proj4string=CRS("+init=epsg:2393"), ...)
      return(invisible(.self))
    },
    
    loadBoundary = function(thin=TRUE, tolerance=0.1) {
      loadBoundaryGADM(country="FIN", thin=thin, tolerance=tolerance)
      return(invisible(.self))
    }
  )
)

RussiaStudyArea <- setRefClass(
  Class = "RussiaStudyArea",
  contains = "StudyArea",
  methods = list(
    initialize = function(...) {
      callSuper(region="Russia", proj4string=CRS("+init=epsg:3410"), ...)
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
      context$getFileName(dir=context$processedDataDirectory, name="Boundary", response=response, region=region)
    },
    
    saveBoundary = function(fileName=getBoundaryFileName(), study) {
      i <- RussiaWTCIntersections$new(study=study)
      i$loadIntersections()
      districts <- loadDistricts()
      boundary <<- subset(districts, NAME_LAT %in% i$intersections$District_Lat & ADM4_ID %in% i$intersections$RegionID)
      save(boundary, file=fileName)
      return(invisible(.self))
    },
    
    loadBoundary = function(thin=TRUE, tolerance=0.1) {
      load(fileName=getBoundaryFileName(), envir=as.environment(.self))
      
      #loadBoundaryGADM(country="RUS", thin=thin, tolerance=tolerance)      
      
      return(invisible(.self))
    },
    
    loadHabitatRaster = function() {
      return(invisible(.self))
    }
  )
)
