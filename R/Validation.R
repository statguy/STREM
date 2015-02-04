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
    
    #getTracksFileIterations = function() {
    #  tracks <- SimulatedTracks(study=study)
    #  tracks$ge
    #},
    
    getIntersectionFileIterations = function() {
      intersections <- SimulatedIntersections(study=study)
      return(intersections$getIntersectionFileIterations())
    },
    
    getPopulationSizeFileIterations = function(modelName) {
      populationSize <- SimulationPopulationSize(study=study, modelName=modelName)
      return(populationSize$getPopulationSizeFileIterations())
    },
    
    getEstimatesFileIterations = function(modelName) {
      model <- SimulatedSmoothModelSpatioTemporal(study=study, modelName=modelName)
      return(model$getEstimatesFileIterations())
    },

    validateTemporalIntersection = function() {
      iterations <- getIntersectionFileIterations()
      if (length(iterations) == 0) return(NULL)
      
      intersections <- ldply(iterations, function(iteration) {
        message("Processing iteration ", iteration, " of ", length(iterations), " iterations for scenario ", study$response, "...")
        intersections <- study$loadIntersections(iteration=iteration, predictDistances=FALSE)
        x <- intersections$intersections@data
        y <- ddply(x, .(year), function(x) data.frame(intersections=sum(x$intersections)))
        y$iteration <- iteration
        return(y)
      })
      
      intersections$Scenario <- study$response
      
      return(intersections)
    },
    
    getValidationTemporalIntersections = function(scenarios=c("A","B","C","D","E","F")) {
      x <- ldply(scenarios, function(scenario) {
        s <- getStudy(scenario=scenario, isTest=F)
        validation <- Validation(study=s, populationSizeOverEstimate=populationSizeOverEstimate)
        x <- validation$validateTemporalIntersection()
        if (is.null(x) || nrow(x) == 0) return(NULL)
        return(x)
      }, .parallel=TRUE)
      
      return(x)
    },
    
    validateTemporalPopulationSize = function(modelName) {
      library(plyr)
      
      iterations <- getPopulationSizeFileIterations(modelName)
      if (length(iterations) == 0) return(NULL)
      
      populationSize <- ldply(iterations, function(iteration) {
        message("Processing iteration ", iteration, " of ", length(iterations), " iterations for scenario ", study$response, "...")
        populationSize <- study$loadPopulationSize(iteration=iteration, modelName=modelName)
        x <- populationSize$sizeData
        
        if (is.null(x) || nrow(x) == 0 || any(x$Estimated > populationSizeOverEstimate)) {
          message("Estimation failed for iteration ", iteration, " scenario ", study$response)
          return(NULL)
        }
        x$iteration <- iteration
        return(x)
      }, .parallel=TRUE)
            
      if (is.null(populationSize)) return(NULL)
      populationSize$Scenario <- study$response
      populationSize$modelName <- modelName
      return(populationSize)
    },
        
    getValidationTemporalPopulationSizes = function(scenarios=c("A","B","C","D","E","F"), modelNames) {
      x <- ddply(expand.grid(scenario=scenarios, modelName=modelNames, stringsAsFactors=FALSE), .(scenario, modelName), function(x) {
        s <- getStudy(scenario=x$scenario, isTest=F)
        validation <- Validation(study=s, populationSizeOverEstimate=populationSizeOverEstimate)
        x <- validation$validateTemporalPopulationSize(x$modelName)
        if (is.null(x) || nrow(x) == 0) return(NULL)
        gc()
        return(x)
      })
      
      return(x)
    },
    
    summarizePopulationSize = function(populationSizes) {
      library(plyr)
      
      ddply(populationSizes, .(Scenario, Year, modelName), function(x, probs) {
        y <- data.frame(n=nrow(x), Estimated=mean(x$Estimated), Observed=mean(x$Observed))
        q <- quantile(x$Estimated, probs=probs)
        y$Estimated.q1 <- q[1]
        y$Estimated.q2 <- q[2]
        q <- quantile(x$Observed, probs=probs)
        y$Observed.q1 <- q[1]
        y$Observed.q2 <- q[2]
        return(y)
      }, probs=c(.0025, .975))
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
      
      iterations <- getEstimatesFileIterations(modelName)
      spatialCorrelation <- data.frame()
      if (length(iterations) == 0) return(spatialCorrelation)
      
      for (iteration in iterations) {
        message("Iteration ", iteration, " of ", length(iterations), " iterations of scenario ", study$response, "...")
        
        estimates <- study$getModel(modelName=modelName, iteration=iteration)
        estimates <- study$loadEstimates(estimates)
        estimates$collectEstimates()
        
        estimated <- estimates$getPopulationDensityInterpolate(templateRaster=template, maskPolygon=NULL, getSD=F)
        if (is.null(estimated$mean)) {
          message("Estimation failed for iteration ", iteration, " scenario ", study$response, ".")
          next
        }
        
        estimated$mean$rasterStack <- stack(mask(estimated$mean$rasterStack, coverArea) * coverArea)
        
        tracks <- study$loadTracks(iteration=iteration)
        cc <- ddply(subset(tracks$tracks, yday == 29), .(year, id), function(x) x[1,])[,c("year","x","y","herdSize")]
        x <- data.frame()  
        
        for (i in unique(tracks$tracks$year) - min(tracks$tracks$year) + 1) {
          year0 <- i + min(tracks$tracks$year) - 1
          
          message("Iteration ", iteration, " of ", length(iterations), " iterations, year ", year0, ", scenario ", study$response, "...")
          
          subcc <- subset(cc, year == year0, select=c("x","y","herdSize"))
          true <- rasterize(subcc[,c("x","y")], template, field=subcc$herdSize, fun='sum', background=0)
          true <- mask(true, coverArea)
          
          correlation <- cor(estimated$mean$rasterStack[[i]][], true[], use="complete.obs", method="spearman")
          estpop <- sum(estimated$mean$rasterStack[[i]][], na.rm=T)
          truepop <- sum(true[], na.rm=T)
          realtrue <- subset(tracks$truePopulationSize, Year == year0)$Observed
                    
          x <- rbind(x, data.frame(Year=year0, Correlation=correlation, True=truepop, Estimated=estpop, TrueError=realtrue-truepop))
        }
        
        if (any(x$Estimated > populationSizeOverEstimate)) {
          message("Estimation failed for iteration ", iteration, " scenario ", study$response, ".")
          next
        }
        
        x$iteration <- iteration
        spatialCorrelation <- rbind(spatialCorrelation, x)
      }
      
      spatialCorrelation$Scenario <- study$response
      spatialCorrelation$modelName <- modelName
      
      return(spatialCorrelation)      
    },
    
    getValidationSpatialPopulationSizes = function(scenarios=c("A","B","C","D","E","F"), modelNames) {
      y <- expand.grid(scenario=scenarios, modelName=modelNames, stringsAsFactors=FALSE)
      spatialCorrelations <- data.frame()
      for (i in 1:nrow(y)) {
        x <- y[i,]
        s <- getStudy(scenario=x$scenario, isTest=F)
        validation <- Validation(study=s, populationSizeOverEstimate=populationSizeOverEstimate)
        z <- validation$validateSpatialPopulationSize(x$modelName)
        if (nrow(z) == 0) next
        spatialCorrelations <- rbind(spatialCorrelations, z)
      }
      return(spatialCorrelations)
      
      #x <- ddply(expand.grid(scenario=scenarios, modelName=modelNames, stringsAsFactors=FALSE), .(scenario, modelName), function(x) {
      #  s <- getStudy(scenario=x$scenario, isTest=F)
      #  validation <- Validation(study=s, populationSizeOverEstimate=populationSizeOverEstimate)
      #  x <- validation$validateSpatialPopulationSize(x$modelName)
      #  gc()
      #  if (nrow(x) == 0) return(NULL)
      #  return(x)
      #}, .parallel=T)
      
      return(x)
    },
    
    getCredibilityIntervalsValidationFileName = function(modelName, iteration) {
      return(study$context$getLongFileName(study$context$scratchDirectory, name="CIValidation", response=study$response, region=study$studyArea$region, tag=paste(modelName, iteration, sep="-")))
    },
    
    validateCredibilityIntervals = function(modelName, iteration, nSamples=100, populationSizeOverEstimate=1e99, save=F) {      
      model <- study$getModel(modelName=modelName, iteration=iteration)
      model$modelName <- modelName
      model$loadEstimates()
      model$collectEstimates()
      
      tracks <- study$loadTracks(iteration=iteration)
      habitatWeights <- study$getHabitatWeights(tracks=tracks, iteration=iteration, save=F)
      
      posteriorSamples <- model$samplePosterior(n=nSamples)
      populationSize <- data.frame()
      
      for (sample in 1:nSamples) {
        message("Processing sample ", sample, " / ", nSamples, " of iteration ", iteration, "...")
        
        model$data <- posteriorSamples[[sample]]
        model$data$fittedMean <- model$data$z
        model$data$year <- model$data$t
        populationDensity <- model$getPopulationDensity(getSD=FALSE)
        
        x <- model$getPopulationSize(populationDensity=populationDensity$mean, tracks=tracks, habitatWeights=habitatWeights)$sizeData
        if (any(x$Estimated >= populationSizeOverEstimate)) {
          message("Estimation failed for iteration ", iteration, ".")
          next
        }
        x$sample <- sample
        populationSize <- rbind(populationSize, x)
      }

      populationSize$scenario <- study$response
      populationSize$modelName <- modelName
      populationSize$iteration <- iteration
      
      if (save) {
        fileName <- getCredibilityIntervalsValidationFileName(modelName=modelName, iteration=iteration)
        message("Saving file to ", fileName, "...")
        save(populationSize, file=fileName)
      }
      
      return(populationSize)  
    },
    
    summarizePopulationSizeCI = function(populationSizeCI, variables=.(scenario, Year), probs=c(.025, .975)) {
      library(plyr)
      
      ddply(populationSizeCI, variables, function(x, probs) {
        y <- data.frame(Estimated=mean(x$Estimated), Observed=mean(x$Observed))
        q <- quantile(x$Estimated, probs=probs)
        y$Estimated.q1 <- q[1]
        y$Estimated.q2 <- q[2]
        return(y)
      }, probs=probs)
    },
    
    getCredibilityIntervalsValidationIterations = function(modelName) {
      return(study$context$getIterationIds(dir=study$context$scratchDirectory, name="CIValidation", response=study$response, region=study$studyArea$region, tag=paste(modelName, "(\\d+)", sep="-")))
    },
    
    loadCredibilityIntervalsValidation = function(modelName, iteration) {
      fileName <- getCredibilityIntervalsValidationFileName(modelName=modelName, iteration=iteration)
      message("Loading file from ", fileName, "...")
      load(fileName)
      return(populationSize)
    },
    
    getValidatedCredibilityIntervalsProportion = function(modelName, probs=c(.025, .975)) {
      library(plyr)
      
      iterations <- getCredibilityIntervalsValidationIterations(modelName=modelName)
      if (length(iterations) == 0) return(data.frame())
      populationSizeCI <- ldply(iterations, function(iteration) {
        loadCredibilityIntervalsValidation(modelName=modelName, iteration=iteration)
      })
      
      x <- summarizePopulationSizeCI(populationSizeCI, variables=.(scenario, Year, iteration), probs=probs)
      validationProportion <- ddply(x, .(scenario), function(x) {
        y <- with(x, Observed >= Estimated.q1 & Observed <= Estimated.q2)
        data.frame(Proportion=mean(y), n=length(y), True=mean(x$Observed), Estimated=mean(x$Estimated), Estimated.q1=mean(x$Estimated.q1), Estimated.q2=mean(x$Estimated.q2))
      })
      validationProportion$modelName <- modelName
      
      return(validationProportion)
    },
    
    getValidatedCredibilityIntervalsProportions = function(scenarios=c("A","B","C","D","E","F"), modelNames, probs=c(.025, .975), probsName="95%") {
      x <- ddply(expand.grid(scenario=scenarios, modelName=modelNames, stringsAsFactors=FALSE), .(scenario, modelName), function(x, probs, probsName) {
        s <- getStudy(scenario=x$scenario, isTest=F)
        validation <- Validation(study=s, populationSizeOverEstimate=populationSizeOverEstimate)
        x <- validation$getValidatedCredibilityIntervalsProportion(x$modelName, probs=probs)
        if (nrow(x) == 0) return(NULL)
        x$Interval <- probsName
        return(x)
      }, probs=probs, probsName=probsName, .parallel=T)
      
      return(x)
    },
    
    # Count also number of tracks
    countFailedEstimations = function(scenarios=c("A","B","C","D","E","F"), modelNames) {
      x <- ddply(expand.grid(scenario=scenarios, modelName=modelNames, stringsAsFactors=FALSE), .(scenario, modelName), function(x) {
        s <- getStudy(scenario=x$scenario, isTest=F)
        validation <- Validation(study=s)
        iterations <- validation$getPopulationSizeFileIterations(x$modelName)
        x <- laply(iterations, function(iteration, modelName, study) {
          message("scenario = ", x$scenario, ", iteration = ", iteration, ", model = ", modelName)
          populationSize <- study$loadPopulationSize(iteration=iteration, modelName=modelName)
          if (any(populationSize$sizeData$Estimated > populationSizeOverEstimate)) return(NA)
          return(iteration)
        }, modelName=x$modelName, study=s, .parallel=T)
        data.frame(nEstimated=length(x[!is.na(x)]))
      })
      return(x)
    }
  )
)
