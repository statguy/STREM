CovariatesContainer <- setRefClass(
  Class = "CovariatesContainer",
  fields = list(
    study = "ANY",
    covariatesId = "character",
    covariateNames = "character"
  ),
  methods = list(
    getCovariatesFileName = function() {
      if (inherits(study, "undefinedField"))
        stop("Parameter 'study' must be specified.")
      if (length(covariatesId) == 0)
        stop("Parameter 'covariatesId' must be specified.")
      return(context$getFileName(dir=study$context$processedDataDirectory, name="Covariates", response=covariatesId, region=study$studyArea$region))
    },
    
    saveCovariates = function(fileName=getCovariatesFileName()) {
      stop("Unimplemented method.")
    },
    
    loadCovariates = function(fileName=getCovariatesFileName()) {
      load(fileName, env=as.environment(.self))
      return(invisible(.self))
    },
    
    setCovariatesId = function(tag) {
      stop("Unimplemented method.")
    },
    
    getSampleLocations = function(...) {
      stop("Unimplemented method.")
    },
    
    associateCovariates = function(...) {
      stop("Unimplemented method.")
    },
    
    addCovariates = function(...) {
      library(plyr)
      covariates <- list(...)
      values <- llply(covariates, function(x) {
        if (!inherits(x, "CovariatesObtainer"))
          stop("Arguments must be of class 'Covariates'.")
        x$preprocess()
        y <- x$extract(getSampleLocations())
        return(y)
      })
      covariateNames <<- do.call(c, lapply(values, colnames))
      do.call(.self$associateCovariates, values)
      return(invisible(.self))
    },
    
    pca = function() {
      
    }
  )
)

###

.interpolateSP = function(xy, z, newXy, transFun=identity, backTransFun=identity) {
  library(sp)
  library(fields)
  if (sum(!is.na(xy[[z]])) == 0) rep(NA, nrow(xy))
  if (sum(!is.na(xy[[z]])) < 10) rep(mean(xy[[z]]), nrow(xy))
  if (proj4string(xy) != proj4string(newXy))
    xy <- sp::spTransform(xy, sp::CRS(sp::proj4string(newXy)))  
  fit <- fields::Tps(sp::coordinates(xy), transFun(xy[[z]]), lambda=0)
  newZ <- as.numeric(predict(fit, sp::coordinates(newXy)))
  return(backTransFun(newZ))
}

# Smoothed single points



.saveRawDiscreteRaster <- function(r, processValues, edgeValues, rowAdd=1000) {
  maskFileName <- tempfile()
  maskFile <- file(maskFileName, "wb")
  edgeFileName <- tempfile()
  edgeFile <- file(edgeFileName, "wb")
  row <- 1
  exit <- FALSE
  while (exit == FALSE) {
    if (row + rowAdd > nrow(r) + 1) {
      rowAdd <- nrow(r) + 1 - row
      exit <- TRUE
    }
    message("Processing row ", row-1, "/", nrow(r), ", reading ", rowAdd, " lines...")
    scanline <- getValues(r, row, rowAdd)
    row <- row + rowAdd
    mask <- scanline %in% processValues
    writeBin(object=mask, con=maskFile, size=1)
    edge <- scanline %in% edgeValues
    writeBin(object=edge, con=edgeFile, size=1)
  }
  close(maskFile)
  close(edgeFile)
  return(list(maskFileName=maskFileName, edgeFileName=edgeFileName))
  #unlink(maskFileName)
  #unlink(edgeFileName)
}

#.smoothDiscreteSubset <- function(r, x, y, kernelFun=.expKernel, scale, processValues, edgeValues) {  

.smoothRaster <- function(r, maskFileName, edgeFileName, kernelFun=.expKernel, scale) {
  rcpp_smoothRaster(maskFileName, edgeFileName, kernelFun, scale, dim(r)[1], dim(r)[2])
}



.identityKernel <- function(size, scale) {
  k <- matrix(0, ncol=2*size+1, nrow=2*size+1)
  k[size+1, size+1] <- 1
  k
}

.expKernel <- function(size, scale) {
  x <- matrix(-size:size, ncol=2*size+1, nrow=2*size+1)
  y <- t(x)
  xy <- sqrt(x^2+y^2)
  xy[xy > size] <- Inf
  k <- exp(-xy / scale)
  k
}

.gaussKernel <- function(size, scale) {
  x <- matrix(-size:size, ncol=2*size+1, nrow=2*size+1)
  y <- t(x)
  xy <- x^2+y^2
  xy[xy > size^2] <- Inf
  k <- exp(-xy / scale)
  k
}

.getKernel <- function(r, scale, kernelFun=.expKernel) {
  resScale <- scale / raster::res(r)[1]
  kernelSize <- round(resScale)
  k <- kernelFun(kernelSize, resScale)
  return(k)  
}

#> x
#coords.x1 
#3693000 
#> y
#coords.x2 
#6979000 

#x <- 3693000
#y <- 6905950
#k <- .getKernel(r=r, scale=100, kernelFun=.expKernel)

.smoothDiscreteSubset <- function(r, x, y, k, scale, processValues, edgeValues) {
  col <- colFromX(r, x)
  row <- rowFromY(r, y)  
  if (is.na(r[row, col]))
    stop("The point is outside the effective area. The smoothing cannot be proceeded.")  
  
  # Cut the kernel if partially outside the effective area
  kernelSize <- (ncol(k)-1)/2 + 1
  startRow <- max(0, row-kernelSize) + 1
  startCol <- max(0, col-kernelSize) + 1
  nrows <- min(dim(r)[1]+1, row+kernelSize) - startRow
  ncols <- min(dim(r)[2]+1, col+kernelSize) - startCol
  xmin <- startRow - (row-kernelSize) - 1
  ymin <- startCol - (col-kernelSize) - 1
  if (!(dim(k)[1] == nrows+xmin & dim(k)[2] == ncols+ymin)) {
    message("Cut kernel ", dim(k)[1], " X ", dim(k)[2], " to ", nrows+xmin, " X ", ncols+ymin)
    k <- k[1:nrows+xmin, 1:ncols+ymin]
  }
  
  # Get edges and process values around the specified point that matches the size of the kernel
  edgeRaster <- raster::getValuesBlock(r, startRow, nrows, startCol, ncols, format='matrix')
  # Set kernel zero at edges for edge correction
  k[edgeRaster %in% edgeValues | is.na(edgeRaster)] <- 0
  # Rescale kernel to constraint smooth values between 0...1
  k <- k / sum(k)
  # Get binary mask of the process (that generates the spatial patterns)
  processRaster <- edgeRaster %in% processValues
  # Find convolution
  smoothValue <- sum(k * processRaster, na.rm=T)
  x <- data.frame(x=x, y=y, scale=scale, value=smoothValue)
  return(x)


  # Get the area of the kernel (could be arbitrarily shaped, non-zero entries indicate the kernel)
  #effectiveArea <- prod(dim(k))
  #effectiveArea <- sum(k > 0)

  # Get edges and process values around the specified point that matches the size of the kernel
  edgeRaster <- raster::getValuesBlock(r, startRow, nrows, startCol, ncols, format='matrix')

  # Set edges to NA
  edgeRaster[edgeRaster %in% edgeValues] <- NA
  # Set area outside the kernel to zero
  edgeRaster <- (k != 0) * edgeRaster
  # Set NA's outside the kernel to zero (edges outside the kernel are ignored too)
  edgeRaster[is.na(edgeRaster) & k == 0] <- 0
  # Count the number of edges (within the area of the kernel)
  edgeCount <- sum(is.na(edgeRaster))
  #edgeCount + sum(edgeRaster>0, na.rm=T) == effectiveArea
  # Get a matrix that indicates the process (that has generated the spatial patterns)
  #processRaster <- (k != 0) * matrix(edgeRaster %in% processValues, nrow=nrow(edgeRaster), ncol=ncol(edgeRaster))
  #processRaster <- (k != 0) * edgeRaster %in% processValues
  processRaster <- matrix(edgeRaster %in% processValues, nrow=nrow(edgeRaster), ncol=ncol(edgeRaster))
#plot(raster(processRaster))
  # Find convolution for the specified point and do edge correction
  smoothValue <- sum(k * processRaster, na.rm=T) * effectiveArea / (effectiveArea - edgeCount)
  
  x <- data.frame(x=x, y=y, scale=scale, value=smoothValue)
  return(x)
}

.smoothDiscreteSubsets <- function(r, coords, kernelFun=.expKernel, scales, processValues, edgeValues, wide=T, .parallel=T) {
  library(raster)
  library(plyr)
  library(sp)
  library(tidyr)
  
  if (missing(r))
    stop("Argument 'r' missing.")
  if (!inherits(r, "RasterLayer"))
    stop("Argument 'r' must be of type 'RasterLayer'")
  if (missing(coords))
    stop("Argument 'coords' missing.")
  
  if (!inherits(r, "RasterLayer"))
    stop("Argument 'r' must of class raster.")
  if (!inherits(coords, "matrix") && !inherits(coords, "data.frame") && !inherits(coords, "SpatialPoints"))
    stop("Argument 'coords' must of class matrix, data.frame or SpatialPoints.")
  
  if (res(r)[1] != res(r)[2])
    stop("Rasters of unequal resolution unsupported.")
  if (nlayers(r) > 1)
    stop("Multiband rasters unsupported. Please supply bands separately.")
  
  if (inherits(coords, "SpatialPoints")) coords <- coordinates(coords)
  
  library(plyr)
  smoothPixels <- plyr::ldply(rev(sort(scales)), function(scale, coords) {
    n.coords <- nrow(coords)
    k <- .getKernel(r=r, scale=scale, kernelFun=.expKernel)
    message("Kernel size = ", dim(k)[1], " X ", dim(k)[2])
    smoothPixels <- ldply(1:n.coords, function(i, coords, n.coords, scale, k) {
    #smoothPixels <- lapply(1:n.coords, function(i, coords, n.coords, scale, k, .parallel) {
      message("Smoothing scale = ", scale, ", for coord = ", i, "/", n.coords, " (", coords[i,1], ",", coords[i,2], ")")
      x <- .smoothDiscreteSubset(r=r, x=coords[i,1], y=coords[i,2], k=k, scale, processValues=processValues, edgeValues=edgeValues)
      return(x)
    }, coords=coords, n.coords=n.coords, scale=scale, k=k, .parallel=.parallel) # Inner loop parallel strategy slower for small kernels, but faster for big kernels
    return(smoothPixels)
  }, coords=coords)
  
  if (wide) smoothPixels <- tidyr::spread(smoothPixels, scale, value)
  return(smoothPixels)
}

CovariatesObtainer <- setRefClass(
  Class = "CovariatesObtainer",
  fields = list(
    study = "Study",
    covariateId = "character"
  ),
  methods = list(
    getCacheFileName = function() {
      if (inherits(study, "uninitializedField"))
        stop("Parameter 'study' must be provided.")
      if (length(covariateId) == 0)
        stop("Parameter 'covariateId' must be provided.")
      return(context$getFileName(dir=study$context$scratchDirectory, name=paste("covariates-cache", covariateId, sep="-"), region=study$studyArea$region))
    },
    
    existCache = function() {
      return(file.exists(getCacheFileName()))
    },
    
    saveCache = function(data) {
      save(data, file=getCacheFileName())
      return(invisible(.self))
    },
    
    loadCache = function() {
      load(getCacheFileName(), env=as.environment(.self))
      return(data)
    },
    
    preprocess = function() {
      stop("Unimplemented abstract method.")
    },
    
    extract = function(xyt) {
      stop("Unimplemented abstract method.")
    }
  )
)

HabitatSmoothCovariates <- setRefClass(
  Class = "HabitatSmoothCovariates",
  contains = "CovariatesObtainer",
  field = list(
    scales = "numeric"
  ),
  methods = list(
    preprocess = function() {
      return(invisible(.self))
    },
    
    extract = function(xyt) {
      library(plyr)
      if (length(scales) == 0)
        stop("Parameter 'scales' must be defined.")
      if (!inherits(xyt, "SpatialPoints"))
        stop("Argument 'xyt' must be of class 'SpatialPoints'.")
      
      # As habitats do not change in time, consider only unique locations
      xy <- as.data.frame(unique(coordinates(xyt)))
      
      habitatWeights <- study$loadHabitatWeights()
      habitatValues <- dlply(habitatWeights$weights[habitatWeights$weights$type != 0,], .(type), function(x) x$habitat)
      names(habitatValues) <- names(habitatWeights$habitatTypes)
      edgeValues <- habitatWeights$weights$habitat[habitatWeights$weights$type == 0]
      
      covariates <- llply(1:length(habitatValues), function(index, edgeValues) {
        message("Processing habitat ", names(habitatValues[index]), "...")
        covariates <- .smoothDiscreteSubsets(r=study$studyArea$habitat, coords=xy, kernelFun=.expKernel, scales=scales, processValues=habitatValues[[index]], edgeValues=edgeValues, .parallel=T)
        covariates <- covariates[,-(1:2),drop=F]
        colnames(covariates) <- paste0("habitat_", names(habitatValues[index]), "_", colnames(covariates))
        return(covariates)
      }, edgeValues=edgeValues)
      
      # Covariates should sum to 1 (approximately).
      # In fact, if all works there is no need to find the last habitat since it will be 1-sum(other habitats). TODO: fix this.
      covariates <- do.call(cbind, covariates)
      sums <- rowSums(covariates, na.rm=T) / length(scales)
      tolerance <- 0.01
      if (any(sums <= 1-tolerance | sums >= 1+tolerance)) {
        print(covariates)
        print(sums)
        stop("Habitat covariates do not sum to 1.")
      }
      
      # Replicate covariates at unique locations over time
      xyz <- cbind(xy, covariates)
      covariates <- join(as.data.frame(coordinates(xyt)), xyz)[,-(1:2),drop=F]
      
      return(covariates)
    }
  )
)

WeatherCovariates <- setRefClass(
  Class = "WeatherCovariates",
  contains = "CovariatesObtainer",
  fields = list(
    apiKey = "character"
  ),
  methods = list(
    initialize = function(...) {
      callSuper(...)
      covariateId <<- "Weather"
    },
    
    preprocess = function() {
      return(invisible(.self))
    },
    
    extract = function(xyt) {
      if (length(apiKey) == 0)
        stop("Parameter 'apiKey' must be specified.")
      if (!inherits(xyt, "SpatialPointsDataFrame"))
        stop("Argument 'xyt' must be of class 'SpatialPointsDataFrame'.")
      
      request <- fmi::FMIWFSRequest$new(apiKey=apiKey)
      allDates <- as.Date(xyt@data[,1])
      uniqueDates <- unique(allDates)
      
      covariates <- list()
      counter <- 1
      for (date in uniqueDates) {
        date <- as.Date(date, origin="1970-1-1")
        message("Finding weather covariates for date ", date, ", ", counter, "/", length(uniqueDates))
        client <- fmi::FMIWFSClient$new(request=request)
        response <- client$getDailyWeather(variables=c("rrday","snow","tday"), startDateTime=date, endDateTime=date, bbox=fmi::getFinlandBBox())
        xytSubset <- xyt[as.Date(xyt$date) == date,]
        rrday <- subset(response, variable == "rrday")
        rrday$measurement[rrday$measurement < 0] <- 0
        rrday <- .interpolateSP(rrday, "measurement", xytSubset, sqrt, square)
        snow <- subset(response, variable == "snow")
        snow$measurement[snow$measurement < 0] <- 0       
        snow <- .interpolateSP(snow, "measurement", xytSubset, sqrt, square)
        tday <- .interpolateSP(subset(response, variable == "tday"), "measurement", xytSubset)
        covariates <- rbind(covariates, data.frame(rrday=rrday, snow=snow, tday=tday))
        counter <- counter + 1
      }

      return(covariates)
    }
  )
)

HumanPopulationDensityCovariates <- setRefClass(
  Class = "HumanPopulationDensityCovariates",
  contains = "CovariatesObtainer",
  fields = list(
    apiKey = "character"
  ),
  methods = list(
    initialize = function(...) {
      callSuper(...)
      covariateId <<- "HumanDensity"
    },
    
    preprocess = function() {
      if (existCache() == FALSE) {
        library(stringr)
        request <- gisfin::GeoStatFiWFSRequest$new()$getPopulationLayers()
        client <- gisfin::GeoStatFiWFSClient$new(request)
        layers <- client$listLayers()
        layers <- stringr::str_subset(layers, "_5km$")
        populationCache <- list()
        for (layer in layers) {
          message("Processing layer ", layer, "...")
          request$getPopulation(layer)
          print(request)
          client <- gisfin::GeoStatFiWFSClient$new(request)
          x <- client$getLayer(layer)
          year <- stringr::str_match(layer, "vaki(\\d+)_")[2]
          y <- sp::spTransform(x, sp::CRS(study$studyArea$proj4string))
          z <- sp::SpatialPixelsDataFrame(points=sp::coordinates(y), data=y@data[,"vaesto",drop=F], proj4string=sp::CRS(y@proj4string), tolerance=0.7)
          populationCache[year] <- z
        }
        saveCache(populationCache)
      }      
      return(invisible(.self))
    },
    
    extract = function(xyt) {
      if (existCache() == FALSE)
        stop("Cache does not exist. Must run preprocess() first.")
      if (!inherits(xyt, "SpatialPointsDataFrame"))
        stop("Argument 'xyt' must be of class 'SpatialPointsDataFrame'.")

      populationCache <- loadCache()
      availYears <- as.numeric(names(populationCache))
      years <- as.POSIXlt(xyt@data[,1])$year + 1900
      uniqueYears <- unique(years)
      covariates <- list()
      counter <- 1
      for (year in uniqueYears) {
        message("Finding population covariates for year ", year, ", ", counter, "/", length(uniqueYears))
        xytSubset <- xyt[years == year,]
        closestYearIndex <- which.min(abs(availYears - year))
        population <- xytSubset %over% populationCache[[closestYearIndex]]
        population[is.na(population)] <- 0
        covariates <- rbind(covariates, data.frame(human=population$vaesto))
        counter <- counter + 1
      }
      
      return(covariates)
    }
  )
)




### TODO: fix
FinlandCovariates <- setRefClass(
  Class = "FinlandCovariates",
  fields = list(
    study = "Study",
    covariatesName = "character",
    covariates = "data.frame"
  ),
  methods = list(
    #initialize = function(study, covariatesName, ...) {
    #  callSuper(study=study, covariatesName=covariatesName, ...)
    #  return(invisible(.self))
    #},
    
    # Note: Currently uses the same raster file for all inheriting classes.
    # Can be changed by adding covariatesName field in the file name.
    getPopulationDensityFileName = function(year) {
      if (inherits(study, "uninitializedField"))
        stop("Provide study parameters.")
      study$context$getFileName(dir=study$context$scratchDirectory, name="PopulationDensityRaster", response=year, region=study$studyArea$region, ext="")
    },
    
    savePopulationDensityYear = function(year, aggregationFactor=1) {
      library(ST)
      library(raster)
      
      x <- getStatFiPopulationDensityRaster(year=year, aggregationFactor=aggregationFactor) # The original grid is 1 x 1 km^2
      names(x) <- paste("year", year, sep="")
      brick(stack(x), filename=getPopulationDensityFileName(year), overwrite=TRUE)
      
      return(invisible(.self))
    },
    
    cachePopulationDensity = function(years, aggregationFactor=1) {      
      for (year in years) {
        message("Caching year ", year, "...")
        savePopulationDensityYear(year=year, aggregationFactor=aggregationFactor)
      }
      return(invisible(.self))
    },
    
    loadPopulationDensityYear = function(year) {
      library(raster)
      return(brick(getPopulationDensityFileName(year)))
    },
    
    getPopulationDensityCovariates = function(xyt) {
      library(plyr)
      library(raster)
      library(ST)

      if (!all(c("id", "year") %in% names(xyt)))
        stop("Missing variables in the data.")
      
      populationDensity <- ddply(raster::as.data.frame(xyt), .(year), function(x, proj4string) {
        year <- x$year[1]
        message("Finding population density covariates for year ", year, "...")
        
        populationDensityRaster <- loadPopulationDensityYear(year=year)
        xy <- SpatialPoints(cbind(x$x, x$y), proj4string=CRS(proj4string))
        populationDensity <- as.vector(getStatFiPopulationDensity(xy, populationDensityRaster))
        
        return(data.frame(id=x$id, populationDensity=populationDensity))
      }, proj4string=proj4string(xyt))
      
      return(populationDensity)
    },
    
    # Note: Currently uses the same raster file for all inheriting classes.
    # Can be changed by adding covariatesName field in the file name.
    getWeatherFileName = function(year) {
      if (inherits(study, "uninitializedField"))
        stop("Provide study parameters.")
      study$context$getFileName(dir=study$context$scratchDirectory, name="WeatherRaster", response=year, region=study$studyArea$region, ext="")
    },
    
    saveWeatherYear = function(year, templateRaster, fmiApiKey) {
      library(ST)
      library(raster)
      library(plyr)
      
      query <- "fmi::observations::weather::daily::multipointcoverage"
      startTime <- paste(year, "-01-01T00:00:00Z", sep="")
      endTime <- paste(year, "-03-31T23:59:59Z", sep="")
      parameters <- list(startTime=startTime, endTime=endTime, bbox="19.0900,59.3000,31.5900,70.130", crs="EPSG::4326")
      weatherStationFile <- queryFMIData(apiKey=fmiApiKey, queryStored=query, parameters=parameters)
      weather <- parseFMIWeatherStationMultipointCoverage(weatherStationFile, study$studyArea$proj4string)
      weather <- as.data.frame(weather)
      weather <- data.frame(weather, breakDownDate(weather$date))
      
      weatherAggregated <- ddply(weather, .(x, y, month), function(x) {
        data.frame(x=x$x[1], y=x$y[1], year=x$year[1], month=x$month[1], rrday=mean(x$rrday, na.rm=T), tday=mean(x$tday, na.rm=T), snow=mean(x$snow, na.rm=T), tmin=mean(x$tmin, na.rm=T), tmax=mean(x$tmax, na.rm=T))
      })
      
      interpolateVariable <- function(weather, variable, templateRaster, transform=identity, inverseTransform=identity) {
        message("Interpolating variable = ", variable, "...")
        weather <- multiRasterInterpolate(xyz=weather[,c("x","y",variable,"month")],
          variables=.(month),
          templateRaster=templateRaster,
          transform=transform,
          inverseTransform=inverseTransform)
        names(weather) <- paste("month", 1:3, variable, sep="")
        return(weather)
      }
      
      x <- list(interpolateVariable(weatherAggregated, "rrday", templateRaster, sqrt, square),
        interpolateVariable(weatherAggregated, "tday", templateRaster),
        interpolateVariable(weatherAggregated, "snow", templateRaster, sqrt, square))
      brick(stack(x), filename=getWeatherFileName(year), overwrite=TRUE)
      
      return(invisible(.self))
    },
    
    loadWeatherYear = function(year) {
      library(raster)
      return(brick(getWeatherFileName(year)))
    },
    
    cacheWeather = function(years, fmiApiKey) {
      for (year in years) 
        saveWeatherYear(year, study$getTemplateRaster(), fmiApiKey=fmiApiKey)
      return(invisible(.self))
    },
    
    getWeatherCovariates = function(xyt) {
      library(raster)
      library(plyr)
      library(sp)
      
      if (!all(c("id", "date") %in% names(xyt)))
        stop("Missing variables in the data.")
      
      x <- data.frame(raster::as.data.frame(xyt), breakDownDate(xyt$date))
      weatherCovariates <- ddply(x, .(year), function(x, proj4string) {
        year <- x$year[1]
        message("Finding weather covariates for year ", year, "...")  
        weather <- loadWeatherYear(year=year)
        
        points <- SpatialPoints(x[,c("x","y")], proj4string=CRS(proj4string))
        weatherValues <- raster::extract(weather, points)
        x$month[x$month > 3] <- 3

        getColumn = function(variable) {
          columns <- paste("month", x$month, variable, sep="")
          index <- cbind(seq_along(columns), match(columns, colnames(weatherValues)))
          return(weatherValues[index])
        }
        
        y <- x[,c("id"),drop=F]
        y$rrday <- getColumn("rrday")
        y$snow <- getColumn("snow")
        y$tday <- getColumn("tday")
        
        return(y)
      }, proj4string=proj4string(xyt), .parallel=TRUE)
      
      return(weatherCovariates)
    },
    
    getCovariatesFileName = function(name) {
      if (inherits(study, "uninitializedField"))
        stop("Provide study parameter.")
      if (length(covariatesName) == 0)
        stop("Provide covariatesName parameter.")
      return(study$context$getFileName(dir=study$context$processedDataDirectory, name=covariatesName, response=study$response, region=study$studyArea$region))
    },
    
    saveCovariates = function(xyt, cache=FALSE, fmiApiKey, impute=FALSE, save=TRUE) {
      library(dplyr)
      
      if (!inherits(xyt, "SpatialPoints"))
        stop("xyt must be a class of SpatialPointsDataFrame")
      
      if (!all(c("id", "year") %in% names(xyt)))
        stop("Missing variables in the data.")
      
      years <- sort(unique(xyt$year))
      if (cache) cachePopulationDensity(years=years, aggregationFactor=4)
      populationDensityCovariates <- getPopulationDensityCovariates(xyt)
      
      if (cache) cacheWeather(years=years, fmiApiKey=fmiApiKey)
      weatherCovariates <- getWeatherCovariates(xyt)
      
      covariates <<- data.frame(populationDensityCovariates, select(weatherCovariates, -c(year,id)))
      #covariates <<- data.frame(populationDensityCovariates, covariates[,!colnames(covariates) %in% c("id","year")])
      if (impute) imputeMissingCovariates(coordinates(xyt))
      
      if (save) save(covariates, file=getCovariatesFileName(covariatesName))

      return(invisible(.self))
    },
    
    loadCovariates = function() {
      load(getCovariatesFileName(), envir=as.environment(.self))
      return(invisible(.self))
    },
    
    imputeMissingCovariates = function(xy) {
      xyt <- data.frame(xy, covariates)
      xyt <- ddply(xyt, .(year), function(x) {
        x <- inverseDistanceWeightningImpute(x, "rrday")
        x <- inverseDistanceWeightningImpute(x, "snow")
        x <- inverseDistanceWeightningImpute(x, "tday")
        return(x)        
      })
      xyt <- arrange(xyt, year, id)
      
      if (sum(xyt$rrday < 0 | xyt$rrday > 250) != 0) {
        message("Invalid rrday:")
        print(xyt[xyt$rrday < 0 | xyt$rrday > 250,])
      }
      if (sum(xyt$snow < 0 | xyt$snow > 220) != 0) {
        message("Invalid snow:")
        print(xyt[xyt$snow < 0 | xyt$snow > 220,])
      }
      if (sum(xyt$tday < -55 | xyt$tday > 40)) {
        message("Invalid tday:")
        print(xyt[xyt$tday < -55 | xyt$tday > 40,])
      }
      
      xyt$x <- NULL
      xyt$y <- NULL
      covariates <<- xyt
      return(invisible(.self))
    }
  )
)

