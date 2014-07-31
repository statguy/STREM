Validation <- setRefClass(
  "Validation",
  fields = list(
    study = "Study",
    populationSizeOverEstimate = "numeric"
  ),
  methods = list(
    initialize = function(..., populationSizeOverEstimate=2000) {
      callSuper(...)
      populationSizeOverEstimate <<- populationSizeOverEstimate
    },
    
    getPopulationSizeFileIterations = function(modelName) {
      populationSize <- SimulationPopulationSize(study=study, modelName=modelName)
      return(populationSize$getPopulationSizeFileIterations())
    },
    
    getEstimatesFileIterations = function(modelName) {
      model <- SimulatedSmoothModel(study=study, modelName=modelName)
      return(model$getEstimatesFileIterations())
    },
    
    getFailedEstimationRate = function(modelName) {
      # TODO
    },
    
    validateTemporalPopulationSize = function(modelName) {
      #iterations <- study$context$getIterationIds(dir=study$context$resultDataDirectory, name="PopulationSize", response=study$response, region=study$studyArea$region, tag=paste("(\\d+)", modelName, sep="-"))
      iterations <- getPopulationSizeFileIterations(modelName)
      size <- data.frame()
      for (iteration in iterations) {
        message("Processing iteration ", iteration, " of ", length(iterations), " iterations.")
        populationSize <- study$loadPopulationSize(iteration=iteration, modelName=modelName)
        populationSize$sizeData$Iteration <- iteration
        if (any(populationSize$sizeData$Estimated > populationSizeOverEstimate)) {
          message("Estimation failed for iteration ", iteration, ".")
          next
        }
        x <- populationSize$sizeData
        x$iteration <- iteration
        size <- rbind(size, x)
      }
      
      size$scenario <- scenario
      size$modelName <- modelName
      
      return(size)
    },
        
    validateSpatialPopulationSize = function(modelName) {
      library(raster)
      
      # 100 km x 100 km grid
      template <- raster(extent(study$studyArea$boundary), nrows=12, ncols=6, crs=study$studyArea$proj4string)
      template@extent@xmax <- template@extent@xmin + dim(template)[2] * 100 * 1e3
      template@extent@ymax <- template@extent@ymin + dim(template)[1] * 100 * 1e3
      coverArea <- rasterize(study$studyArea$boundary, template, getCover=T)
      coverArea[coverArea[] == 0] <- NA
      coverArea <- coverArea / 100
      
      #iterations <- study$context$getIterationIds(dir=study$context$scratchDirectory, name=modelName, response=study$response, region=study$studyArea$region, tag="(\\d+)")
      iterations <- getEstimatesFileIterations(modelName)
      spatialCorrelation <- data.frame()
      
      for (iteration in iterations) {
        message("Iteration ", iteration, " of ", length(iterations), " iterations.")
        
        estimates <- SimulatedSmoothModelSpatioTemporal(study=study, iteration=iteration)
        estimates$modelName <- modelName
        estimates <- study$loadEstimates(estimates)
        estimates$collectEstimates()
        estimated <- estimates$getPopulationDensity(templateRaster=template, maskPolygon=NULL, getSD=F)
        estimated$mean$rasterStack <- stack(mask(estimated$mean$rasterStack, coverArea) * coverArea)
        
        tracks <- study$loadTracks(iteration=iteration)
        cc <- ddply(subset(tracks$tracks, yday == 29), .(year, id), function(x) x[1,])[,c("year","x","y")]
        x <- data.frame()  
        
        for (i in unique(tracks$tracks$year) - min(tracks$tracks$year) + 1) {
          year0 <- i + min(tracks$tracks$year) - 1
          true <- rasterize(subset(cc, year == year0, select=c("x","y")), template, field=1, fun='count', background=0)
          true <- mask(true, coverArea)
          
          correlation <- cor(estimated$mean$rasterStack[[i]][], true[], use="complete.obs", method="spearman")
          estpop <- sum(estimated$mean$rasterStack[[i]][], na.rm=T)
          truepop <- sum(true[], na.rm=T)
          
          x <- rbind(x, data.frame(Year=year0, Correlation=correlation, True=truepop, Estimated=estpop))
        }
        
        if (any(x$Estimated) > populationSizeOverEstimate) {
          message("Estimation failed for iteration ", iteration, ".")
          next
        }
        
        x$iteration <- iteration
        spatialCorrelation <- rbind(spatialCorrelation, x)
      }
      
      spatialCorrelation$scenario <- scenario
      spatialCorrelation$modelName <- modelName
      
      return(spatialCorrelation)      
    },
    
    getCredibilityIntervalsValidationFileName = function(modelName, iteration) {
      return(study$context$getLongFileName(study$context$scratchDirectory, name=modelName, response=study$response, region=study$studyArea$region, tag=iteration))
    },
    
    validateCredibilityIntervals = function(modelName, iteration, nSamples=100, save=F) {      
      model <- SimulatedSmoothModelSpatioTemporal$new(study=study, iteration=iteration)
      model$modelName <- modelName
      model$loadEstimates()
      model$collectEstimates()
      
      tracks <- study$loadTracks(iteration=iteration)
      # TODO: habitat weights
      
      posteriorSamples <- model$samplePosterior(n=nSamples)
      populationSize <- data.frame()
      
      for (sample in 1:nSamples) {
        message("Processing sample ", sample, " / ", nSamples, " of iteration ", iteration, "...")
        
        populationDensity <- model$getPopulationDensity.internal(posteriorSamples[[sample]])
        x <- model$getPopulationSize(populationDensity=populationDensity, tracks=tracks)$sizeData
        if (any(x$Estimated > populationSizeOverEstimate)) {
          message("Estimation failed for iteration ", iteration, ".")
          break
        }
        x$sample <- sample
        populationSize <- rbind(populationSize, x)
      }

      populationSize$scenario <- scenario
      populationSize$modelName <- modelName
      populationSize$iteration <- iteration
      
      if (save) {
        save(populationSize, file=getCredibilityIntervalsValidationFileName(modelName=modelName, iteration=iteration))
      }
      
      return(populationSize)  
    },
    
    loadCredibilityIntervalsValidation = function(modelName, iteration) {
      load(getCredibilityIntervalsValidationFileName(modelName=modelName, iteration=iteration))
      return(populationSize)
    }
  )
)
