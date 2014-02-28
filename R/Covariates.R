FinlandCovariates <- setRefClass(
  Class = "FinlandCovariates",
  fields = list(
    study = "Study",
    populationDensityCache = "RasterStack",
    covariatesName = "character",
    covariates = "data.frame"
  ),
  methods = list(
    initialize = function(...) {
      callSuper(...)
      return(invisible(.self))
    },
    
    cachePopulationDensity = function(years, aggregationFactor=1) {
      library(ST)
      library(raster)
      
      for (year in years) {
        message("Caching year ", year, "...")
        populationDensity <- getStatFiPopulationDensityRaster(year=year, aggregationFactor=aggregationFactor) # The original grid is 1 x 1 km^2
        names(populationDensity) <- paste("year", year, sep="")
        populationDensityCache <<- if (inherits(populationDensityCache, "uninitializedField")) stack(populationDensity)
        else stack(populationDensityCache, populationDensity)
      }
    },
    
    getPopulationDensityCovariates = function(xyt) {
      library(plyr)
      library(raster)
      library(ST)

      if (!c("year","id","x","y") %in% names(xyt))
        stop("Missing variables in the data.")
      
      populationDensity <- ddply(as.data.frame(xyt), .(year), function(x, proj4string) {
        year <- x$year[1]
        message("Processing year ", year, "...")
        xy <- SpatialPoints(cbind(x$x, x$y), proj4string=CRS(proj4string))
        
        index <- paste("year", year, sep="")
        if (!index %in% names(populationDensityCache))
          stop("No corresponding raster for year ", year)

        populationDensityRaster <- populationDensityCache[[index]]
        populationDensity <- getStatFiPopulationDensity(xy, populationDensityRaster)
        return(data.frame(id=x$id, populationDensity=populationDensity))
      }, proj4string=proj4string(xyt), .parallel=TRUE)
      
      return(populationDensity)
    },
    
    getWeatherFileName = function(year) {
      if (inherits(study, "uninitializedField"))
        stop("Provide study parameters.")
      study$context$getFileName(dir=study$context$scratchDirectory, name="WeatherRaster", response=year, region=study$studyArea$region, ext="")
    },
    
    saveWeatherYear = function(year, templateRaster, apiKey) {
      library(ST)
      library(raster)
      
      query <- "fmi::observations::weather::daily::multipointcoverage"
      startTime <- paste(year, "-01-01T00:00:00Z", sep="")
      endTime <- paste(year, "-03-31T23:59:59Z", sep="")
      parameters <- list(startTime=startTime, endTime=endTime, bbox="19.0900,59.3000,31.5900,70.130", crs="EPSG::4326")
      weatherStationFile <- queryFMIData(apiKey=apiKey, queryStored=query, parameters=parameters)
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
    
    cacheWeather = function(years, apiKey) {
      library(foreach)
      foreach(year=years) %dopar% saveWeatherYear(year, study$getTemplateRaster(), apiKey=apiKey)
      return(invisible(.self))
    },
    
    getWeatherCovariates = function(xyt) {
      library(raster)
      library(plyr)
      
      if (!c("id","date","x","y") %in% names(xyt))
        stop("Missing variables in the data.")
      
      x <- data.frame(as.data.frame(xyt), breakDownDate(xyt$date))
      weatherCovariates <- ddply(x, .(year), function(x, proj4string) {
        year <- x$year[1]
        message("Processing year ", year, "...")        
        weather <- loadWeatherYear(year=year)
        
        points <- SpatialPoints(x[,c("x","y")], proj4string=CRS(proj4string))
        weatherValues <- extract(weather, points)
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
      }, proj4string=proj4string(xyt), .parallel=TRUE, .inform=T)
      
      return(weatherCovariates)
    },
    
    getCovariatesFileName = function(name) {
      if (inherits(study, "uninitializedField"))
        stop("Provide study parameters.")
      return(study$context$getFileName(dir=study$context$resultDataDirectory, name=covariatesName, response=study$response, region=study$studyArea$region))
    },
    
    saveCovariates = function(xyt, apiKey) {
      library(raster)
      
      if (!c("id","year","x","y") %in% names(xyt))
        stop("Missing variables in the data.")
      years <- sort(unique(xyt$year))
      
      cacheWeather(years=years, apiKey=apiKey)
      weatherCovariates <- getWeatherCovariates(xyt, aggregationFactor=4)
      
      cachePopulationDensity(years=years)
      populationDensityCovariates <- getPopulationDensityCovariates(xyt)
      
      covariates <<- merge(populationDensityCovariates, weatherCovariates, sort=FALSE)
      save(covariates, file=getCovariatesFileName(covariatesName))
      
      populationDensityCache <<- stack()
      gc()
      
      return(invisible(.self))
    },
    
    loadCovariates = function() {
      load(getCovariatesFileName(), envir=as.environment(.self))
      return(invisible(.self))
    }
  )
)

