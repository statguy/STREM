library(INLA)
setOldClass("inla.mesh")
setOldClass("inla.spde2")
setOldClass("inla.data.stack")
setOldClass("inla")

Model <- setRefClass(
  Class = "Model",
  fields = list(
    study = "Study",
    data = "data.frame",
    locations = "matrix",
    #covariates = "data.frame",
    model = "formula",
    coordsScale = "numeric",
    years = "integer",
    mesh = "inla.mesh",
    spde = "inla.spde2",
    index = "list",
    A = "Matrix",
    obsStack = "inla.data.stack",
    predStack = "inla.data.stack",
    fullStack = "inla.data.stack",
    result = "inla",
    node = "list",
    modelName = "character",
    offsetScale = "numeric"
  ),
  methods = list(
    initialize = function(...) {
      callSuper(...)
      return(invisible(.self))
    },
    
    getEstimatesFileName = function() {
      if (inherits(study, "undefinedField") | length(modelName) == 0)
        stop("Provide study and modelName parameters.")
      return(study$context$getFileName(study$context$resultDataDirectory, name=modelName, response=study$response, region=study$studyArea$region))
    },
    
    saveEstimates = function(fileName=getEstimatesFileName()) {
      message("Saving result to ", fileName, "...")
      save(locations, data, model, coordsScale, years, mesh, spde, index, A, fullStack, result, offsetScale, file=fileName)
    },
    
    loadEstimates = function(fileName=getEstimatesFileName()) {
      load(fileName, envir=as.environment(.self))
      return(invisible(.self))
    },
    
    estimate = function(save=FALSE, fileName=getEstimatesFileName()) {
      stop("Unimplemented method.")
    },
    
    getUnscaledMesh = function() {
      mesh.unscaled <- mesh
      mesh.unscaled$loc <- mesh.unscaled$loc / coordsScale
      return(mesh.unscaled)
    },
    
    getUnscaledMeshCoordinates = function() {
      return(getUnscaledMesh()$loc[,1:2])
    }
  )
)

SmoothModel <- setRefClass(
  Class = "SmoothModel",
  contains = "Model",
  fields = list(
    family = "character",
    interceptPrior = "list",
    rhoPrior = "list"
  ),
  methods = list(
    initialize = function(...) {
      callSuper(modelName="SmoothModel", ...)
      return(invisible(.self))
    },
    
    setupInterceptPrior = function(priorParams) {
      # TODO: Find area from raster ignoring NA values if habitat weights are used (raster area sligtly different from the boundary area)
      # TODO: Link function is fixed, change if needed
      
      area <- study$studyArea$boundary@polygons[[1]]@area
      toTheta <- function(x, area) log(x / area)
      fromTheta <- function(x, area) exp(x) * area
      
      populationSizeMean <- priorParams$mean
      populationSizeDeviation <- priorParams$sd
      
      x <- toTheta(populationSizeMean, area)
      y <- transformDeviation(populationSizeMean, populationSizeDeviation, toTheta, area=area)
      interceptPrior <<- list(mean=x, prec=1/y)
      message("Intercept prior mean = ", interceptPrior$mean, ", sd = ", 1/interceptPrior$prec)
      
      z <- fromTheta(qnorm(c(0.025,0.5,0.975), x, y), area)
      message("Intercept prior 0.025 0.5 0.975 quantiles = ", signif(z[1],2), " ", signif(z[2],2), " ", signif(z[3],2))
      
      return(invisible(.self))
    },

    setupRhoPrior = function(priorParams) {
      warning("Prior for rho not in use.")
      rhoPrior <<- list(theta=priorParams)
      #z <- qnorm(c(0.025, 0.5, 0.975), priorParams$param[1], priorParams$param[2])
      #message("Rho prior 0.025 0.5 0.975 quantiles = ", signif(z[1],2), " ", signif(z[2],2), " ", signif(z[3],2))
      return(invisible(.self))
    },
    
    getObservedOffset = function(distance=data$distance) {
      return(2/pi * data$length * data$duration * distance / offsetScale)
    },
    
    getPredictedOffset = function(distance=mean(data$distance)) {
      # Put on approximately the same scale with observed crossings
      return(rep(2/pi * 12000 * 1 * distance, mesh$n * length(years) / offsetScale))
    },
    
    setup = function(intersections, meshParams, offsetScale=1, family="nbinomial") {
      library(INLA)
      library(plyr)
      
      offsetScale <<- offsetScale
      data <<- intersections$getData()
      
      p <- ddply(data, .(year), function(x) {
        data.frame(pzero=sum(x$intersections==0)/length(x$intersections), mean=mean(x$intersections))
      })
      message("Proportion of zeros, mean each year:")
      print(p)
      message("Proportion of intersections:")
      print(prop.table(table(data$intersections)))
      
      message("Constructing mesh...")
      
      coordsScale <<- meshParams$coordsScale
      locations <<- intersections$getCoordinates() * coordsScale
      mesh <<- inla.mesh.create.helper(points.domain=locations,
                                       min.angle=meshParams$minAngle,
                                       max.edge=meshParams$maxEdge * coordsScale,
                                       cutoff=meshParams$cutOff * coordsScale)

      years <<- as.integer(sort(unique(data$year)))
      nYears <- length(years)
      groupYears <- as.integer(data$year - min(years) + 1)
            
      spde <<- inla.spde2.matern(mesh)
      index <<- inla.spde.make.index("st", n.spde=spde$n.spde, n.group=nYears)
      
      family <<- family
      model <<- response ~ -1 + intercept + f(st, model=spde, group=st.group, control.group=list(model="ar1"))#, hyper=rhoPrior))
      A <<- inla.spde.make.A(mesh, loc=locations, group=groupYears, n.group=nYears)
      
      message("Number of nodes in the mesh = ", mesh$n)
      message("Average survey route length = ", mean(data$length))
      message("Average count duration = ", mean(data$duration))
      message("Average animal track length = ", mean(data$distance))
      
      obsStack <<- inla.stack(data=list(response=data$intersections,
                                        E=getObservedOffset(),
                                        link=1),
                               A=list(A),
                               effects=list(c(index, list(intercept=1))),
                               tag="observed")
      
      predStack <<- inla.stack(data=list(response=NA,
                                         E=1, #getPredictedOffset(),
                                         link=1),
                               A=list(1),
                               effects=list(c(index, list(intercept=1))),
                               tag="predicted")
            
      return(invisible(.self))
    },
    
    estimate = function(save=FALSE, fileName=getEstimatesFileName(), verbose=TRUE) {
      library(INLA)
      
      fullStack <<- inla.stack(obsStack, predStack)
      #stackData <- inla.stack.data(fullStack)
      stackData <- inla.stack.data(fullStack, spde=spde)
      
      control.fixed <- if (!inherits(interceptPrior, "undefinedField"))
         list(mean=list(intercept=interceptPrior$mean), prec=list(intercept=interceptPrior$prec))
      else NULL
      
      message("Estimating population density...")
      
      result <<- inla(model,
                      family=family,
                      #data=inla.stack.data(fullStack, spde=spde),
                      data=stackData,
                      E=stackData$E,
                      verbose=verbose,
                      control.fixed=control.fixed,
                      control.predictor=list(A=inla.stack.A(fullStack), link=stackData$link, compute=TRUE),
                      control.compute=list(cpo=FALSE, dic=TRUE))
      
      if (is.null(result$ok) | result$ok == FALSE) {
        warning("INLA failed to run.")
      }
      else {
        if (save) saveEstimates(fileName=fileName)
      }
    },

    collectEstimates = function(observationWeights=1, predictionWeights=1, predictAtNodesOnOriginalScale=FALSE) {
      library(INLA)
      
      message("Processing fitted values...")
      
      indexObserved <- inla.stack.index(fullStack, "observed")$data
      data$eta <<- result$summary.linear.predictor$mean[indexObserved] + log(observationWeights)
      data$fittedMean <<- result$summary.fitted.values$mean[indexObserved] * observationWeights
      data$fittedSD <<- result$summary.fitted.values$sd[indexObserved] * observationWeights^2
      
      #stackData <- inla.stack.data(fullStack, tag="observed")
      #observedOffset <- stackData$E[indexObserved]
      observedOffset <- getObservedOffset()
      
      message("Fitted values sums all years:")
      message("observed = ", sum(data$intersections))
      message("estimated = ", sum(data$fittedMean * observedOffset))
      
      message("Processing random effects...")
      
      indexPredicted <- inla.stack.index(fullStack, "predicted")$data
      node$mean <<- matrix(result$summary.fitted.values$mean[indexPredicted] * predictionWeights, nrow=mesh$n, ncol=length(years))
      node$sd <<- matrix(result$summary.fitted.values$sd[indexPredicted] * predictionWeights^2, nrow=mesh$n, ncol=length(years))
      node$spatialMean <<- matrix(result$summary.random$st$mean, nrow=mesh$n, ncol=length(years))
      node$spatialSd <<- matrix(result$summary.random$st$sd, nrow=mesh$n, ncol=length(years))
      #predictedOffset <- matrix(stackData$E[indexPredicted], nrow=mesh$n, ncol=length(years)) * matrix(getPredictedOffset(), nrow=mesh$n, ncol=length(years))
      if (predictAtNodesOnOriginalScale) predictedOffset <- matrix(getPredictedOffset(), nrow=mesh$n, ncol=length(years))
      
      stat <- data.frame()
      for (year in years) {
        yearWhich <- data$year == year
        yearIndex <- year - min(years) + 1
        x <- if (predictAtNodesOnOriginalScale) 
          data.frame(
            Year=years[yearIndex],
            
            Observed=sum(data$intersections[yearWhich]),
            EstimatedAtObserved=sum(data$fittedMean[yearWhich] * observedOffset[yearWhich]),
            EstimatedAtNodes=sum(node$mean[,yearIndex] * predictedOffset[,yearIndex]),
            
            ObservedScaled=sum(data$intersections[yearWhich] / observedOffset[yearWhich]),
            EstimatedAtObservedScaled=sum(data$fittedMean[yearWhich]),
            EstimatedAtNodesScaled=sum(node$mean[,yearIndex]),
            EstimatedAtNodesScaledAdjusted=sum(node$mean[,yearIndex]) * nrow(locations[yearWhich,]) / mesh$n,
            
            ObservedOffset=mean(observedOffset[yearWhich]),
            predictedOffset=mean(predictedOffset[,yearIndex])
          )
        else
          data.frame(
            Year=years[yearIndex],
            
            Observed=sum(data$intersections[yearWhich]),
            EstimatedAtObserved=sum(data$fittedMean[yearWhich] * observedOffset[yearWhich]),
            
            ObservedScaled=sum(data$intersections[yearWhich] / observedOffset[yearWhich]),
            EstimatedAtObservedScaled=sum(data$fittedMean[yearWhich]),
            EstimatedAtNodesScaled=sum(node$mean[,yearIndex]),
            EstimatedAtNodesScaledAdjusted=sum(node$mean[,yearIndex]) * nrow(locations[yearWhich,]) / mesh$n,
            
            ObservedOffset=mean(observedOffset[yearWhich])
          )
        
        stat <- rbind(stat, x)
      }
      
      message("Year by year summary:")
      print(stat)
      message("Column sums:")
      print(colSums(stat))
      message("Column means:")
      print(colMeans(stat))
      message("Correlations:")
      print(cor(stat[!colnames(stat) %in% c("Year")]))
      
      return(invisible(stat))
    },
    
    collectHyperparameters = function() {
      library(INLA)
      
      # TODO: Handle offset scale
      warning("TODO: Handle offset scale")
      
      message("Processing hyperparameters...")
      spdeResult <- inla.spde2.result(result, "st", spde)
      logKappa <- getINLAEstimates(spdeResult$marginals.log.kappa[[1]], coordsScale=1)
      logTau <- getINLAEstimates(spdeResult$marginals.log.tau[[1]], coordsScale=1) ## scale ??
      sd <- getINLAEstimates(spdeResult$marginals.log.variance.nominal[[1]], fun=function(x) sqrt(exp(x)))
      range <- getINLAEstimates(spdeResult$marginals.range.nominal[[1]], coordsScale=coordsScale)
      # range = sqrt(8)/kappa
      # sd = 1/(4*pi*kappa^2*tau^2)
      
      x <- rbind(logKappa=logKappa,
                 logTau=logTau,
                 sd=sd,
                 range=range)
      if (any(rownames(result$summary.hyperpar)=="GroupRho for st"))
        x <- rbind(x, rho=result$summary.hyperpar["GroupRho for st",])
      
      return(x)
    },
    
    # TODO: fix?
    DEPRECATED_getPredictedIntersections = function() {
      x <- data
      indexObserved <- inla.stack.index(fullStack, "observed")$data
      x$offset <- inla.stack.data(fullStack)$E[indexObserved]
      intersections <- ddply(x, .(year), function(x) {
        data.frame(Observed=sum(x$intersections), Predicted=sum(x$fittedMean * x$offset), logObsOff=log(sum(x$intersections / x$offset)), eta=sum(x$eta))
      })
      return(intersections)
    },
    
    project = function(projectValues, projectionRaster=study$getTemplateRaster(), maskPolygon, weights=1) {
      library(INLA)
      library(raster)
      
      projector <- inla.mesh.projector(getUnscaledMesh(),
                                       dims=c(ncol(projectionRaster), nrow(projectionRaster)),
                                       xlim=c(xmin(projectionRaster), xmax(projectionRaster)),
                                       ylim=c(ymin(projectionRaster), ymax(projectionRaster)))
      projectedEstimates <- inla.mesh.project(projector, projectValues)
      values(projectionRaster) <- t(projectedEstimates[,ncol(projectedEstimates):1]) * weights
      
      if (!missing(maskPolygon))
        projectionRaster <- mask(projectionRaster, maskPolygon)
      
      return(invisible(projectionRaster))
    },
    
    getPopulationDensity = function(templateRaster=study$getTemplateRaster(), maskPolygon=study$studyArea$boundary, getSD=TRUE) {
      if (length(node) == 0)
        stop("Did you forgot to run collectEstimates() first?")
      
      library(raster)

      meanPopulationDensityRaster <- SpatioTemporalRaster$new(study=study)
      sdPopulationDensityRaster <- SpatioTemporalRaster$new(study=study)
      
      xyzMean <- data.frame(Year=data$year, locations / coordsScale, z=data$fittedMean)
      xyzSD <- data.frame(Year=data$year, locations / coordsScale, z=data$fittedSD)
      cellArea <- prod(res(templateRaster)) # m^2
      
      for (year in sort(unique(xyzMean$Year))) {
        yearIndex <- year - min(xyzMean$Year) + 1
        message("Processing year ", year, "...")
        
        meanRaster <- project(projectValues=node$mean[,yearIndex] / offsetScale, projectionRaster=templateRaster, maskPolygon=maskPolygon, weights=cellArea)
        meanPopulationDensityRaster$addLayer(meanRaster, year)
        
        if (getSD) {
          sdRaster <- project(projectValues=node$sd[,yearIndex] / offsetScale^2, projectionRaster=templateRaster, maskPolygon=maskPolygon, weights=cellArea^2)
          sdPopulationDensityRaster$addLayer(sdRaster, year)
        }
      }
      
      return(invisible(list(mean=meanPopulationDensityRaster, sd=sdPopulationDensityRaster)))
    },
    
    associateMeshLocationsWithDate = function() {
      library(plyr)
      library(fields)
      
      #if (!inherits(xyt, "SpatialPointsDataFrame"))
      #  stop("Argument xyt must be of class SpatialPointsDataFrame.")
      #if (any(!c("year","date") %in% names(xyt)))
      #  stop("Argument xyt must have data columns names year, date.")
      #if (inherits(mesh, "uninitializedField"))
      #  stop("Run setupMesh() first.")
      
      nYears <- length(years)
      meshLocations <- getUnscaledMesh()$loc[,1:2]
      nMeshLocations <- nrow(meshLocations)
      #predictLocations <- repeatMatrix(meshLocations, nYears)      
      #id <- rep(1:nMeshLocations, nYears)
      #year <- rep(years, each=nMeshLocations)
      
      xyt <- ddply(data, .(year), function(x) {
        message("Finding the closest observations for mesh nodes for year ", x$year[1], "...")
        xyt <- data.frame()
        for (i in 1:nMeshLocations) {
          xy <- meshLocations[i,,drop=F]
          d <- rdist(xy, x[,1:2] / coordsScale)
          xyt <- rbind(xyt, data.frame(x=xy[,1], y=xy[,2], date=x[which.min(d),"date"], id=i))
        }
        return(xyt)
      })
      
      xyt <- SpatialPointsDataFrame(xyt[,c("x","y")],
                                    data.frame(id=xyt$id, year=xyt$year, date=xyt$date),
                                    proj4string=study$studyArea$proj4string)
      return(xyt)
    },
    
    saveMeshNodeCovariates = function(save=FALSE) {
      stop("Method saveMeshNodeCovariates unimplemented.")
    },
    
    plotMesh = function(surveyRoutes) {
      library(INLA)
      library(sp)
      meshUnscaled <- mesh
      meshUnscaled$loc <- mesh$loc / coordsScale
      plot(meshUnscaled, asp=1, main="", col="gray")
      plot(study$studyArea$boundary, border="black", add=T)
      
      if (missing(surveyRoutes)) {}
      # TODO: plot intersection coordinates
      else plot(surveyRoutes$surveyRoutes, col="blue", add=T)
      
      return(invisible(.self))
    },
    
    plotProjection = function(projectionRaster=study$getTemplateRaster(), variable, weights=1) {
      meanRaster <- project(projectValues=variable, projectionRaster=projectionRaster, weights=weights)
      op <- par(mar=rep(0, 4))
      plot(meanRaster)
      plot(study$studyArea$boundary, add=T)
      points(locations / coordsScale, col=data$intersections)
      par(op)
    },
    
    plotSpatialEstimatedMean = function(yearIndex=1) {
      plotProjection(variable=node$mean[,yearIndex], weights=prod(res(study$getTemplateRaster())))
    },
    
    plotSpatialEstimatedSD = function(yearIndex=1) {
      plotProjection(variable=node$sd[,yearIndex], weights=prod(res(study$getTemplateRaster()))^2)
    },

    plotSpatialRandomEffectMean = function(yearIndex=1) {
      plotProjection(variable=node$spatialMean[,yearIndex])
    },
    
    plotSpatialRandomEffectSD = function(yearIndex=1) {
      plotProjection(variable=node$spatialSd[,yearIndex])
    },
    
    plotTemporal = function(observationWeights=getObservedOffset()) {
      library(ggplot2)
      library(reshape2)
      
      indexObserved <- inla.stack.index(fullStack, "observed")$data
      stackData <- inla.stack.data(fullStack, tag="observed")
      observedOffset <- stackData$E[indexObserved] * observationWeights
      
      stat <- data.frame()
      for (year in years) {
        yearWhich <- data$year == year
        yearIndex <- year - min(years) + 1
        x <- data.frame(
          Year=years[yearIndex],
          Observed=sum(data$intersections[yearWhich]),
          Estimated=sum(data$fittedMean[yearWhich] * observedOffset[yearWhich])
        )
        stat <- rbind(x, stat)
      }
      
      x <- melt(stat, id="Year")
      p <- ggplot(x, aes(Year, value, group=variable, color=variable)) + geom_line()
      print(p)
      
      return(invisible(p))
    }
  )
)

SimulatedSmoothModel <- setRefClass(
  Class = "SimulatedSmoothModel",
  contains = c("SmoothModel"),
  fields = list(
    iteration = "integer"
  ),
  methods = list(
    getEstimatesFileName = function() {
      if (inherits(study, "undefinedField") | length(modelName) == 0 | length(iteration) == 0)
        stop("Provide study, modelName and iteration parameters.")
      return(study$context$getLongFileName(study$context$scratchDirectory, name=modelName, response=study$response, region=study$studyArea$region, tag=iteration))
    },
    
    getPopulationSize = function(tracks, withHabitatWeights=FALSE) {
      if (length(node) == 0)
        stop("Run collectEstimates() first.")
      populationDensity <- getPopulationDensity(getSD=FALSE)
      
      if (withHabitatWeights) {
        habitatWeights <- CORINEHabitatWeights$new(study=study)
        if (missing(tracks)) tracks <- study$loadTracks(iteration=iteration)
        habitatSelection <- tracks$getHabitatPreferences(habitatWeightsTemplate=habitatWeights, nSamples=30, save=FALSE) # TODO: save
        habitatWeights$setHabitatSelectionWeights(habitatSelection)
        habitatWeighstRaster <- habitatWeights$getWeightsRaster(save=FALSE)
        populationDensity$mean$weight(habitatWeightsRaster)
      }
      
      populationSize <- populationDensity$mean$integrate(volume=SimulationPopulationSize$new(study=study, iteration=iteration))
      if (missing(tracks)) populationSize$loadValidationData
      else populationSize$loadValidationData(tracks=tracks)
      
      return(invisible(populationSize))
    }
  )
)

SimulatedSmoothModelNoOffset <- setRefClass(
  Class = "SimulatedSmoothModelNoOffset",
  contains = c("SimulatedSmoothModel"),
  fields = list(
    iteration = "integer"
  ),
  methods = list(
    getObservedOffset = function(distance) {
      return(rep(1, nrow(data)))
    },
    
    getPredictedOffset = function(distance) {
      return(rep(1, mesh$n * length(years)))
    }
  )
)

FinlandSmoothModel <- setRefClass(
  Class = "FinlandSmoothModel",
  contains = c("SmoothModel", "FinlandCovariates"),
  fields = list(
  ),
  methods = list(
    initialize = function(...) {
      callSuper(modelName="SmoothModel", covariatesName="FinlandSmoothModelCovariates", ...)
      return(invisible(.self))
    },
    
    getPredictedOffset = function(distance=covariates$distance) {
      if (length(distance) == 0) {
        saveMeshNodeCovariates()
        distance <- covariates$distance
        #tracks <- study$loadTracks()
        #distance <- tracks$getMeanDistance()
      }
      #if (length(distance) == 1) distance <- rep(distance, mesh$n * length(years))
      return(2/pi * 12000 * 1 * distance * offsetScale)
    },
    
    saveMeshNodeCovariates = function(save=FALSE) {
      meshNodeLocationsSP <- associateMeshLocationsWithDate()
      saveCovariates(meshNodeLocationsSP, impute=TRUE, save=save)
      covariates$distance <<- predictDistances()
      return(invisible(.self))
    },
    
    predictDistances = function(formula=study$getDistanceCovariatesModel(), intervalH=study$getTrackSampleInterval()) {
      distances <- study$predictDistances(formula=formula, data=covariates, intervalH=intervalH)
      return(distances)
    }
  )
)

FinlandSmoothModelNoDistances <- setRefClass(
  Class = "FinlandSmoothModelNoDistances",
  contains = c("SmoothModel", "FinlandCovariates"),
  fields = list(
  ),
  methods = list(
    initialize = function(...) {
      callSuper(modelName="SmoothModel", covariatesName="FinlandSmoothModelCovariates", ...)
      return(invisible(.self))
    },

    getObservedOffset = function(distance) {
      return(rep(1, nrow(data)))
    },
    
    getPredictedOffset = function(distance) {
      return(rep(1, mesh$n * length(years)))
    }
  )
)

FinlandRussiaSmoothModel <- setRefClass(
  Class = "FinlandRussiaSmoothModel",
  contains = c("SmoothModel"),
  fields = list(
  ),
  methods = list(
    initialize = function(...) {
      callSuper(modelName="SmoothModel", ...)
      return(invisible(.self))
    }
  )
)
