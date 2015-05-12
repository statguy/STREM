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
    }  
  )
)

###

.interpolateSP = function(xy, z, newXy, transFun=identity, backTransFun=identity) {
  library(sp)
  library(fields)
  if (proj4string(xy) != proj4string(newXy))
    xy <- sp::spTransform(xy, sp::CRS(sp::proj4string(newXy)))
  fit <- fields::Tps(sp::coordinates(xy), transFun(xy[[z]]), lambda=0)
  newZ <- predict(fit, sp::coordinates(newXy))
  return(backTransFun(newZ))
}

# Smoothed single points

.identityKernel <- function(size, scale) {
  k <- matrix(0, ncol=2*size+1, nrow=2*size+1)
  k[size+1, size+1] <- 1
  k
}

.expKernel <- function(size, scale) {
  x <- matrix(-size:size, ncol=2*size+1, nrow=2*size+1)
  y <- t(x)
  k <- exp(-sqrt(x^2+y^2) / scale)
  k
}

.gaussKernel <- function(size, scale) {
  x <- matrix(-size:size, ncol=2*size+1, nrow=2*size+1)
  y <- t(x)
  k <- exp(-(x^2+y^2) / scale)
  k
}

.smoothSubset <- function(r, x, y, kernelFun=.expKernel, scale) {
  col <- colFromX(r, x)
  row <- rowFromY(r, y)  
  if (is.na(r[row, col]))
    stop("The point is outside the effective area. The smoothing cannot be proceeded.")  
  
  # Construct the full kernel
  resScale <- scale / res(r)[1]
  kernelSize <- round(resScale)
  k <- kernelFun(kernelSize, resScale)
  fullArea <- prod(dim(k))
  kernelSize <- kernelSize + 1
  
  # Cut the kernel if partially outside the effective area
  startRow <- max(0, row-kernelSize) + 1
  startCol <- max(0, col-kernelSize) + 1
  nrows <- min(dim(r)[1]+1, row+kernelSize) - startRow
  ncols <- min(dim(r)[2]+1, col+kernelSize) - startCol
  xmin <- startRow - (row-kernelSize) - 1
  ymin <- startCol - (col-kernelSize) - 1
  k <- k[1:nrows+xmin, 1:ncols+ymin]
  k <- k / sum(k) # Scale the kernel to produce convoluted values between 0...1
  effectiveArea <- prod(dim(k))
  
  # Get the process values around the point that matches the size of the kernel
  rasterSubset <- getValuesBlock(r, startRow, nrows, startCol, ncols, format='matrix')
  # Count the number of edges
  edgeCount <- sum(is.na(rasterSubset))
  # Smooth and do edge correction
  smoothValue <- sum(k * rasterSubset, na.rm=T) * effectiveArea / (effectiveArea - edgeCount)
  
  return(data.frame(x=x, y=y, scale=scale, value=smoothValue))
}

.smoothSubsets <- function(r, coords, kernelFun=expKernel, scales, .parallel=T) {
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
  
  library(plyr)
  smoothPixels <- ldply(scales, function(scale, coords) {
    n.coords <- if (inherits(coords, "SpatialPoints")) length(coords)
    else nrow(coords)
    smoothPixels <- ldply(1:n.coords, function(i, coords, n.coords, scale) {
      message("Smoothing scale = ", scale, ", for coord = ", i, "/", n.coords, " (", coords[i,1], ",", coords[i,2], ")")
      x <- .smoothSubset(r=r, x=coords[i,1], y=coords[i,2], kernelFun=kernelFun, scale=scale)
      return(x)
    }, coords=coords, n.coords=n.coords, scale=scale)
    return(smoothPixels)
  }, coords=coords, .parallel=.parallel)
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
      load(getCacheFileName())
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
        rrday <- .interpolateSP(subset(response, variable == "rrday"), "measurement", xytSubset, sqrt, square)
        snow <- .interpolateSP(subset(response, variable == "snow"), "measurement", xytSubset, sqrt, square)
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

SmoothHabitatCovariates <- setRefClass(
  Class = "SmoothHabitatCovariates",
  contains = "CovariatesObtainer",
  field = list(
  ),
  methods = list(
    extract = function(xyt) {
      transects <- study$loadSurveyRoutes(findLengths=F)
      habitatWeights <- study$loadHabitatWeights()
      habitatValues <- dlply(habitatWeights$weights[habitatWeights$weights$type != 0,], .(type), function(x) x$habitat)
      edgeValues <- habitatWeights$weights$habitat[habitatWeights$weights$type == 0]
      covariates <- smoothSubsets(study$studyArea$habitat, coords=transects$centroids, scales=scales, processValues=habitatValues, edgeValues=edgeValues)
      if (save) saveCovariates()
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

