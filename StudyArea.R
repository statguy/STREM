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
    coordinateScale = "integer"
  ),
  methods = list(
    initialize = function(context=NA, region=NA, proj4string=NA) {
      callSuper(context=context, region=region, proj4string=proj4string)      
    },
    
    prepareHabitatRaster = function(rawRasterFile) {
      if (missing(rawRasterFile))
        stop("rawRasterFile or boundary parameter missing.")
      saveHabitatRasterFile(rawRasterFile=rawRasterFile)      
    },
    
    maskLargeRaster = function(inputRasterFile, outputRasterFile) {
      require(maptools)
      
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
    },
    
    getHabitatRasterFile = function() {
      stop("Override this method.")
    },
    
    loadHabitatRaster = function() {
      if (!file.exists(getHabitatRasterFile()))
        stop("Habitat raster file ", getHabitatRasterFile(), " not found.")
      habitat <<- raster(getHabitatRasterFile())
      projection(habitat) <<- attr(proj4string, "projargs")
    },
            
    findLargestPolygon = function(gadm) {
      require(maptools)
      p <- list()
      for (i in 1:length(gadm)) {
        largest <- which.max(sapply(gadm@polygons[[i]]@Polygons, slot, "area"))
        p[[i]] <- Polygons(list(Polygon(gadm@polygons[[i]]@Polygons[[largest]])), ID=i)
      }
      p <- SpatialPolygonsDataFrame(SpatialPolygons(p, proj4string=gadm@proj4string), data=gadm@data)
      return(p)
    },
    
    thinGADM = function(gadm, tolerance) {
      p <- findLargestPolygon(gadm)
      return(thinnedSpatialPoly(p, tolerance=tolerance))
    },    
    
    loadBoundary = function(country, level=0, subregion, mainland=FALSE, thin=FALSE, coordinateScale=as.integer(1), tolerance=5000) {  
      coordinateScale <<- coordinateScale
      
      require(rgdal)
      
      if (thin==TRUE & mainland==TRUE) mainland <- FALSE
      
      gadm <- raster::getData('GADM', country=country, level=level, path=context$scratchDirectory)
      
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
    },
    
    getHabitatFrequenciesFileName = function() {
      return(file.path(context$processedDataDirectory, paste("HabitatFrequencies-", region, ".RData", sep="")))
    },
    
    saveHabitatFrequencies = function(habitatRaster) {
      habitatFrequencies <- freq(habitatRaster, progress="text")
      save(habitatFrequencies, file=getHabitatFrequenciesFileName())
    },
    
    loadHabitatFrequencies = function() {
      if (!file.exists(getHabitatFrequenciesFileName()))
        stop("Habitat frequencies file ", getHabitatFrequenciesFileName(), " not found.")
      load(getHabitatFrequenciesFileName())
      return(habitatFrequencies)
    },
    
    plotHabitatRaster = function() {
      plot(habitat)
      plot(boundary, add=T)
    }
  )
)

TestStudyArea <- setRefClass(
  Class = "TestStudyArea",
  contains = "StudyArea",
  methods = list(
    initialize = function(...) {
      callSuper(region="test", proj4string=CRS("+init=epsg:2393"), ...)
            
      x <- -5
      habitat <<- SpatialPolygons(list(Polygons(list(Polygon(matrix(c(325+x,700, 325+x,712, 338+x,712, 338+x,700, 325+x,700)*10000, ncol=2, byrow=T))), ID=1)), proj4string=proj4string)
      
      loadHabitatRaster()
      habitat <<- crop(habitat, boundary)
      habitat <<- readAll(habitat)
      habitat[is.na(habitat[])] <<- 44
    },
    
    getHabitatRasterFile = function() {
      return(file.path(context$scratchDirectory, paste("HabitatRaster-Finland.tif")))
    }
  )
)

FinlandStudyArea <- setRefClass(
  Class = "FinlandStudyArea",
  contains = "StudyArea",
  methods = list(
    initialize = function(...) {
      callSuper(region="Finland", proj4string=CRS("+init=epsg:2393"), ...);
    },
    
    loadBoundary = function(country="FIN", thin=TRUE, tolerance=0.1) {
      callSuper(country=country, thin=thin, tolerance=tolerance)
    },
    
    getHabitatRasterFile = function() {
      return(file.path(context$scratchDirectory, paste("HabitatRaster-Finland.tif")))
    }
  )
)
