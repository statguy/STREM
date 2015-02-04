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
    model = "formula",
    coordsScale = "numeric",
    years = "integer",
    result = "inla",
    modelName = "character",
    offsetScale = "numeric",
    family = "character"
  ),
  methods = list(
    initialize = function(...) {
      callSuper(...)
      return(invisible(.self))
    },
    
    setModelName = function(family, randomEffect) {
      modelName <<- paste("SmoothModel", family, randomEffect, sep="-")
      return(invisible(.self))
    },
    
    getEstimatesFileName = function() {
      if (inherits(study, "undefinedField") | length(modelName) == 0)
        stop("Provide study and modelName parameters.")
      return(study$context$getFileName(study$context$resultDataDirectory, name=modelName, response=study$response, region=study$studyArea$region))
    },
    
    saveEstimates = function(fileName=getEstimatesFileName()) {
      stop("Unimplemented method saveEstimates.")
    },
    
    loadEstimates = function(fileName=getEstimatesFileName()) {
      message("Loading estimates from ", fileName)
      load(fileName, envir=as.environment(.self))
      return(invisible(.self))
    },
        
    getUnscaledMesh = function() {
      mesh.unscaled <- mesh
      mesh.unscaled$loc <- mesh.unscaled$loc / coordsScale
      return(mesh.unscaled)
    },
    
    getUnscaledMeshCoordinates = function() {
      return(getUnscaledMesh()$loc[,1:2])
    },
    
    getUnscaledObservationCoordinates = function() {
      return(locations / coordsScale)
    },
    
    getObservedOffset = function(distance=data$distance) {      
      return(2/pi * data$length * data$duration * distance / offsetScale)
    },
    
    setup = function(intersections, params) {
      library(INLA)
      library(plyr)
      
      if (missing(params))
        stop("Missing params argument.")
      if (!hasMember(params, "family"))
        stop("Missing family parameter.")
      if (!hasMember(params, "timeModel"))
        stop("Missing timeModel parameter.")
      coordsScale <<- if (!hasMember(params, "coordsScale")) 1 else params$coordScale
      offsetScale <<- if (!hasMember(params, "offsetScale")) 1000^2 else params$offsetScale
      family <<- if (!hasMember(params, "family")) "nbinomial" else params$family
      
      setModelName(family=family, randomEffect=params$timeModel)
      data <<- intersections$getData()
      data$response <<- data$intersections
      locations <<- intersections$getCoordinates() * coordsScale
      years <<- as.integer(sort(unique(data$year)))
      model <<- if (hasMember(params, "model")) params$model
      else response ~ 1 + f(year, model=params$timeModel)
      
      return(invisible(.self))      
    },
    
    estimate = function(save=FALSE, fileName=getEstimatesFileName()) {
      stop("Unimplemented method.")
    },
    
    samplePosterior = function(n, index) {
      library(INLA)
      library(plyr)
      posteriorSamples <- inla.posterior.sample(n=n, result=result)
      xy <- getUnscaledObservationCoordinates()
      if (missing(index)) index <- 1:nrow(xy)
      x <- llply(posteriorSamples,
                 function(x, index) {
                   predictorIndex <- grep("Predictor\\.", rownames(x$latent)[index])
                   #data.frame(x=xy[,1], y=xy[,2], z=exp(x$latent[predictorIndex]) / offsetScale, t=data$year)
                   data.frame(x=xy[,1], y=xy[,2], z=exp(x$latent[predictorIndex]), t=data$year)
                 },
                 index=index, .parallel=F)
      return(x)
    },
           
    getPopulationDensity = function(templateRaster=study$getTemplateRaster(), maskPolygon=study$studyArea$boundary, getSD=FALSE) {
      if (is.null(data$fittedMean))
        stop("Did you forgot to run collectEstimates() first?")
      
      library(raster)
      library(plyr)
      
      #area <- maskPolygon@polygons[[1]]@area
      meanPopulationDensity <- ddply(data, .(year), function(x)
        data.frame(z=mean(x$fittedMean / offsetScale), year=x$year[1]))
      
      cellArea <- prod(res(templateRaster))
      meanPopulationDensityRaster <- SpatioTemporalRaster(study=study)$fill(z=meanPopulationDensity$z, boundary=maskPolygon, layerNames=meanPopulationDensity$year, weights=cellArea)
      
      sdPopulationDensityRaster <- if (getSD) {
        sdPopulationDensity <- ddply(data, .(year), function(x)
          data.frame(z=mean(x$fittedSD / offsetScale), year=x$year[1]))
        SpatioTemporalRaster(study=study)$fill(z=meanPopulationDensity$z, boundary=maskPolygon, layerNames=meanPopulationDensity$year, weights=cellArea)
      }
      else SpatioTemporalRaster(study=study)
      
      return(invisible(list(mean=meanPopulationDensityRaster, sd=sdPopulationDensityRaster)))
    },
    
    # Remove this and use *PopulationSize class instead
    getPopulationSize = function(populationDensity, tracks, habitatWeights) {
      if (missing(populationDensity)) populationDensity <- getPopulationDensity(getSD=FALSE)$mean
      if (!missing(habitatWeights)) if (!is.null(habitatWeights)) populationDensity$weight(habitatWeights)
      
      populationSize <- populationDensity$integrate(volume=SimulationPopulationSize(study=study, iteration=iteration, modelName=modelName))
      if (missing(tracks)) populationSize$loadValidationData()
      else populationSize$loadValidationData(tracks=tracks)
      
     return(invisible(populationSize))
    }    
  )
)

AggregatedModel <- setRefClass(
  Class = "AggregatedModel",
  contains = "Model",
  fields = list(
  ),
  methods = list(
    aggregate = function() {
      library(plyr)
      data <<- ddply(data, .(year), function(x) {
        data.frame(response=sum(x$response), intersections=sum(x$intersections), duration=1, length=sum(x$duration*x$length*x$distance), distance=1)
      })      
    },
        
    getPopulationSize = function(populationDensity, tracks, habitatWeights) {
      if (is.null(data$fittedMean))
        stop("Did you forgot to run collectEstimates() first?")
      if (length(unique(data$year)) != nrow(data))
        stop("The data must be aggregated to determine population size")
      area <- if (missing(habitatWeights) || is.null(habitatWeights)) study$studyArea$boundary@polygons[[1]]@area
      else cellStats(habitatWeights, sum) * prod(res(habitatWeights))
      populationSize <- SimulationPopulationSize(study=study, modelName=modelName, iteration=iteration)
      
      for (y in sort(data$year)) {
        size <- subset(data, year == y)$fittedMean * area
        populationSize$addYearSize(y, size)
      }        
      
      if (missing(tracks)) populationSize$loadValidationData()
      else populationSize$loadValidationData(tracks=tracks)
      
      return(populationSize)
    }
  )
)

FMPModel <- setRefClass(
  Class = "FMPModel",
  contains = "AggregatedModel",
  fields = list(
  ),
  methods = list(
    setup = function(intersections, params) {
      coordsScale <<- 1
      offsetScale <<- 1000^2
      modelName <<- "FMPModel"
      
      data <<- intersections$getData()
      data$response <<- data$intersections
      locations <<- intersections$getCoordinates() * coordsScale
      years <<- as.integer(sort(unique(data$year)))
      
      .self$aggregate()
      
      return(invisible(.self))      
    },
    
    saveEstimates = function(fileName=getEstimatesFileName()) {
      message("Saving result to ", fileName, "...")
      save(locations, data, coordsScale, years, offsetScale, file=fileName)
    },
    
    estimate = function(save=FALSE, fileName=getEstimatesFileName(), verbose=TRUE) {
      if (save) saveEstimates(fileName=fileName)
    },
    
    collectEstimates = function(observationWeights=1, predictionWeights=1) {
      index <- 1:nrow(data)

      data$fittedMean <<- data$intersections / getObservedOffset() * observationWeights
      data$fittedSD <<- NA
      observedOffset <- getObservedOffset()
      
      message("Fitted values sums all years:")
      message("observed = ", sum(data$intersections))
      message("estimated = ", sum(data$fittedMean * observedOffset))
      
      stat <- data.frame()
      for (year in years) {
        yearWhich <- data$year == year
        yearIndex <- year - min(years) + 1
        x <- data.frame(
          Year=years[yearIndex],
          Observed=sum(data$intersections[yearWhich]),
          EstimatedAtObserved=sum(data$fittedMean[yearWhich] * observedOffset[yearWhich]),
          ObservedScaled=sum(data$intersections[yearWhich] / observedOffset[yearWhich]),
          EstimatedAtObservedScaled=sum(data$fittedMean[yearWhich]),
          ObservedOffset=mean(observedOffset[yearWhich])
        )
        stat <- rbind(stat, x)
      }
      message("Year by year summary:")
      print(stat)
    }
  )
)

SmoothModelTemporal <- setRefClass(
  Class = "SmoothModelTemporal",
  contains = "Model",
  fields = list(
  ),
  methods = list(
    saveEstimates = function(fileName=getEstimatesFileName()) {
      message("Saving result to ", fileName, "...")
      save(locations, data, model, coordsScale, years, result, offsetScale, family, file=fileName)
    },
            
    estimate = function(save=FALSE, fileName=getEstimatesFileName(), verbose=TRUE) {
      library(INLA)
      
      message("Estimating population density...")
      
      result <<- inla(model,
                      family=family,
                      data=data,
                      E=getObservedOffset(),
                      verbose=verbose,
                      control.predictor=list(compute=TRUE),
                      control.compute=list(cpo=FALSE, dic=TRUE, config=TRUE))
      
      if (is.null(result$ok) | result$ok == FALSE) {
        warning("INLA failed to run.")
      }
      else {
        if (save) saveEstimates(fileName=fileName)
      }
    },
    
    collectEstimates = function(observationWeights=1, predictionWeights=1) {
      library(INLA)
      
      index <- 1:nrow(data)
      data$eta <<- result$summary.linear.predictor$mean[index] + log(observationWeights)
      data$fittedMean <<- result$summary.fitted.values$mean[index] * observationWeights
      data$fittedSD <<- result$summary.fitted.values$sd[index] * observationWeights
      observedOffset <- getObservedOffset()
      
      message("Fitted values sums all years:")
      message("observed = ", sum(data$intersections))
      message("estimated = ", sum(data$fittedMean * observedOffset))      
      
      stat <- data.frame()
      for (year in years) {
        yearWhich <- data$year == year
        yearIndex <- year - min(years) + 1
        x <- data.frame(
          Year=years[yearIndex],
          Observed=sum(data$intersections[yearWhich]),
          EstimatedAtObserved=sum(data$fittedMean[yearWhich] * observedOffset[yearWhich] / offsetScale),
          ObservedScaled=sum(data$intersections[yearWhich] / observedOffset[yearWhich] / offsetScale),
          EstimatedAtObservedScaled=sum(data$fittedMean[yearWhich]),
          ObservedOffset=mean(observedOffset[yearWhich])
        )
        stat <- rbind(stat, x)
      }
      
      message("Year by year summary:")
      print(stat)
      message("Column sums:")
      x <- stat[!colnames(stat) %in% c("Year")]
      print(colSums(x))
      message("Column means:")
      print(colMeans(x))
      message("Correlations:")
      print(cor(x))
      
      return(invisible(stat))
    },
    
    collectHyperparameters = function() {
      library(INLA)
      
      message("Processing hyperparameters...")
      
      x <- data.frame()
      if (any(rownames(result$summary.hyperpar)=="Precision for year")) # RW1, RW2, AR1, ARp, seasonal
        x <- rbind(x, random_effect_precision=result$summary.hyperpar["Precision for year",])
      if (any(rownames(result$summary.hyperpar)=="Rho for year")) # AR1
        x <- rbind(x, rho=result$summary.hyperpar["Rho for year",])
      if (any(rownames(result$summary.hyperpar)=="PACF1 for year")) # ARp
        x <- rbind(x, rho=result$summary.hyperpar["PACF1 for year",])
      if (any(rownames(result$summary.hyperpar)=="PACF2 for year")) # ARp
        x <- rbind(x, rho=result$summary.hyperpar["PACF2 for year",])
      if (any(rownames(result$summary.hyperpar)=="PACF3 for year")) # ARp
        x <- rbind(x, rho=result$summary.hyperpar["PACF3 for year",])
      if (any(rownames(result$summary.hyperpar)=="PACF4 for year")) # ARp
        x <- rbind(x, rho=result$summary.hyperpar["PACF4 for year",])
      
      return(x)
    }
  )
)

SmoothModelMeanTemporal <- setRefClass(
  Class = "SmoothModelMeanTemporal",
  contains = c("SmoothModelTemporal", "AggregatedModel"),
  fields = list(
  ),
  methods = list(
    setModelName = function(family, randomEffect) {
      modelName <<- paste("SmoothModelMean", family, randomEffect, sep="-")
      return(invisible(.self))
    },
    
    setup = function(intersections, params) {
      callSuper(intersections, params)
      .self$aggregate()
      return(invisible(.self)) 
    }
  )
)

SmoothModelSpatioTemporal <- setRefClass(
  Class = "SmoothModelSpatioTemporal",
  contains = "Model",
  fields = list(
    mesh = "inla.mesh",
    spde = "inla.spde2",
    index = "list",
    A = "Matrix",
    augA = "Matrix",
    obsStack = "inla.data.stack",
    augStack = "inla.data.stack",
    predStack = "inla.data.stack",
    fullStack = "inla.data.stack",
    interceptPrior = "list",
    rhoPrior = "list",
    spdeParams = "list",
    node = "list"
  ),
  methods = list(
    saveEstimates = function(fileName=getEstimatesFileName()) {
      message("Saving result to ", fileName, "...")
      save(locations, data, model, coordsScale, years, mesh, spde, index, A, fullStack, result, offsetScale, family, interceptPrior, rhoPrior, spdeParams, file=fileName)
    },
    
    setupInterceptPrior = function(priorParams) {
      # TODO: Find area from raster ignoring NA values if habitat weights are used (raster area sligtly different from the boundary area)
      # TODO: Link function is fixed, change if needed
      
      area <- study$studyArea$boundary@polygons[[1]]@area
      toTheta <- function(x, area) log(x / area)
      fromTheta <- function(x, area) exp(x) * area
      
      populationSizeMean <- priorParams$mean
      populationSizeDeviation <- priorParams$sd
      
      prior.mean <- toTheta(populationSizeMean, area)
      prior.sd <- transformDeviation(populationSizeMean, populationSizeDeviation, toTheta, area=area)
      interceptPrior <<- list(mean=prior.mean, prec=1/prior.sd)

      message("Intercept prior mean = ", prior.mean, ", sd = ", 1/prior.sd)
      
      z <- fromTheta(qnorm(c(0.025,0.5,0.975), prior.mean, prior.sd), area)
      message("Intercept prior 0.025 0.5 0.975 quantiles = ", signif(z[1],2), " ", signif(z[2],2), " ", signif(z[3],2))
      
      return(invisible(.self))
    },
    
    setupTemporalPrior = function(priorParams) {
      rhoPrior <<- list(theta=list(param=c(priorParams$mean, priorParams$sd), initial=priorParams$initial))
      z <- qnorm(c(0.025, 0.5, 0.975), priorParams$mean, priorParams$sd)
      message("Rho prior 0.025 0.5 0.975 quantiles = ", signif(z[1],2), " ", signif(z[2],2), " ", signif(z[3],2))
      return(invisible(.self))
    },
    
    setupSpatialPrior = function(priorParams) {
      kappa.mean <- rangeToKappa(priorParams$range$mean)
      tau.mean <- sigmaToTau(priorParams$sigma$mean, kappa.mean)
      tau.log.sd <- sort(log(tau.mean) - log(sigmaToTau(c(priorParams$sigma$lower, priorParams$sigma$upper), kappa.mean))) / qnorm(c(0.025,0.975))
      kappa.log.sd <- sort(log(kappa.mean) - log(rangeToKappa(c(priorParams$range$lower, priorParams$range$upper)))) / qnorm(c(0.025,0.975))
      
      message("tau mean = ", signif(tau.mean,2), ", sd candidates = ", signif(exp(tau.log.sd[1]),2), " ", signif(exp(tau.log.sd[2]),2), ", sd mean = ", signif(mean(exp(tau.log.sd)),2))
      message("kappa mean = ", signif(kappa.mean,2), ", sd candidates = ", signif(exp(kappa.log.sd[1]),2), " ", signif(exp(kappa.log.sd[2]),2), ", sd mean = ", signif(mean(exp(kappa.log.sd)),2))
      message("log tau mean = ", signif(log(tau.mean),2), ", sd candidates = ", signif(tau.log.sd[1],2), " ", signif(tau.log.sd[2],2), ", sd mean = ", signif(mean(tau.log.sd),2))
      message("log kappa mean = ", signif(log(kappa.mean),2), ", sd candidates = ", signif(kappa.log.sd[1],2), " ", signif(kappa.log.sd[2],2), ", sd mean = ", signif(mean(kappa.log.sd),2))
      message("log tau precision = ", signif(1/mean(tau.log.sd),4), ", log kappa precision = ", signif(1/mean(kappa.log.sd),4))
      
      spdeParams <<- list(kappa=kappa.mean, tau=tau.mean, B.tau=cbind(log(tau.mean), -1, 1), B.kappa=cbind(log(kappa.mean), 0, -1),
                          theta.prior.mean=c(0,0), theta.prior.prec=c(1/mean(tau.log.sd),1/mean(kappa.log.sd)),
                          constr=priorParams$constr)
      return(invisible(.self))
    },
    
    getPredictedOffset = function(distance=mean(data$distance)) {
      # Put on approximately the same scale with observed crossings
      return(rep(2/pi * 12000 * 1 * distance, mesh$n * length(years)) / offsetScale)
    },
    
    setup = function(intersections, params) {
      library(INLA)
      library(R.utils)
      
      callSuper(intersections=intersections, params=params)
      setModelName(family=params$family, randomEffect=paste("matern", params$timeModel, sep="-"))

      if (!hasMember(params, "meshParams"))
        stop("Missing meshParams parameter.")
      
      if (hasMember(params, "interceptPrior")) setupInterceptPrior(params$interceptPrior)
      if (hasMember(params, "temporalPrior")) setupTemporalPrior(params$temporalPrior)
      if (hasMember(params, "spatialPrior")) setupSpatialPrior(params$spatialPrior)
      
      message("Constructing mesh...")
      
      coordsScale <<- params$meshParams$coordsScale
      locations <<- intersections$getCoordinates() * coordsScale
      
      mesh <<- evalWithTimeout(inla.mesh.create.helper(points.domain=locations,
                                                       min.angle=params$meshParams$minAngle,
                                                       max.edge=params$meshParams$maxEdge * coordsScale,
                                                       cutoff=params$meshParams$cutOff * coordsScale),
                               timeout=10, onTimeout="error")
      #mesh <<- inla.mesh.create.helper(points.domain=locations,
      #                                 min.angle=params$meshParams$minAngle,
      #                                 max.edge=params$meshParams$maxEdge * coordsScale,
      #                                 cutoff=params$meshParams$cutOff * coordsScale)

      nYears <- length(years)
      groupYears <- as.integer(data$year - min(years) + 1)
      
      spde <<- if (length(spdeParams) == 0) inla.spde2.matern(mesh)
      else inla.spde2.matern(mesh, B.tau=spdeParams$B.tau, B.kappa=spdeParams$B.kappa, theta.prior.mean=spdeParams$theta.prior.mean, theta.prior.prec=spdeParams$theta.prior.prec, constr=spdeParams$constr)
      index <<- inla.spde.make.index("st", n.spde=spde$n.spde, n.group=nYears)
      
      model <<- if (length(rhoPrior) == 0) response ~ -1 + intercept + f(st, model=spde, group=st.group, control.group=list(model="ar1"))
      else response ~ -1 + intercept + f(st, model=spde, group=st.group, control.group=list(model="ar1", hyper=rhoPrior))
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
      
      boundarySamples <- study$studyArea$sampleBoundary()
      if (!is.null(boundarySamples)) {
        coords <- repeatMatrix(coordinates(boundarySamples), length(years)) * coordsScale
        yearsData <- rep(years, each=length(boundarySamples))      
        groupYears <- as.integer(yearsData - min(yearsData) + 1)
        augA <<- inla.spde.make.A(mesh, loc=coords, group=groupYears, n.group=nYears)
        augStack <<- inla.stack(data=list(response=NA,
                                           E=1,
                                           link=1),
                                 A=list(augA),
                                 effects=list(c(index, list(intercept=1))),
                                 tag="augmented")
      }
      
      predStack <<- inla.stack(data=list(response=NA,
                                         E=1,
                                         link=1),
                               A=list(1),
                               effects=list(c(index, list(intercept=1))),
                               tag="predicted")
      
      return(invisible(.self))
    },
    
    estimate = function(save=FALSE, fileName=getEstimatesFileName(), verbose=TRUE) {
      library(INLA)
      
      fullStack <<- if (is.null(augStack)) inla.stack(obsStack, predStack)
      else inla.stack(obsStack, augStack, predStack)
      stackData <- inla.stack.data(fullStack, spde=spde)
      
      control.fixed <- if (length(interceptPrior) > 0)
         list(mean.intercept=interceptPrior$mean, prec.intercept=interceptPrior$prec)
      else list()
      
      message("Estimating population density...")
      
      result <<- inla(model,
                      family=family,
                      data=stackData,
                      E=stackData$E,
                      verbose=verbose,
                      control.fixed=control.fixed,
                      control.predictor=list(A=inla.stack.A(fullStack), link=stackData$link, compute=TRUE),
                      control.compute=list(cpo=FALSE, dic=TRUE, config=TRUE))
      
      if (is.null(result$ok) | result$ok == FALSE) {
        warning("INLA failed to run.")
      }
      else {
        if (save) saveEstimates(fileName=fileName)
      }
    },

    collectEstimates = function(observationWeights=1, predictionWeights=1, predictAtNodesOnOriginalScale=FALSE) {
      library(INLA)
      
      indexObserved <- inla.stack.index(fullStack, "observed")$data
      data$eta <<- result$summary.linear.predictor$mean[indexObserved] + log(observationWeights)
      data$fittedMean <<- result$summary.fitted.values$mean[indexObserved] * observationWeights
      data$fittedSD <<- result$summary.fitted.values$sd[indexObserved] * observationWeights
      observedOffset <- getObservedOffset()
      
      message("Fitted values sums all years:")
      message("observed = ", sum(data$intersections))
      message("estimated = ", sum(data$fittedMean * observedOffset))
      
      indexPredicted <- inla.stack.index(fullStack, "predicted")$data
      node$mean <<- inla.vector2matrix(result$summary.fitted.values$mean[indexPredicted] * predictionWeights, nrow=mesh$n, ncol=length(years))
      node$sd <<- inla.vector2matrix(result$summary.fitted.values$sd[indexPredicted] * predictionWeights, nrow=mesh$n, ncol=length(years))
      node$spatialMean <<- inla.vector2matrix(result$summary.random$st$mean, nrow=mesh$n, ncol=length(years))
      node$spatialSd <<- inla.vector2matrix(result$summary.random$st$sd, nrow=mesh$n, ncol=length(years))
      if (predictAtNodesOnOriginalScale)
        predictedOffset <- inla.vector2matrix(getPredictedOffset(), nrow=mesh$n, ncol=length(years))
      
      stat <- data.frame()
      for (year in years) {
        yearWhich <- data$year == year
        yearIndex <- year - min(years) + 1
        x <- if (predictAtNodesOnOriginalScale) 
          data.frame(
            Year=years[yearIndex],
            
            Observed=sum(data$intersections[yearWhich]),
            EstimatedAtObserved=sum(data$fittedMean[yearWhich] * observedOffset[yearWhich] / offsetScale),
            EstimatedAtNodes=sum(node$mean[,yearIndex] * predictedOffset[,yearIndex]),
            
            ObservedScaled=sum(data$intersections[yearWhich] / observedOffset[yearWhich] / offsetScale),
            EstimatedAtObservedScaled=sum(data$fittedMean[yearWhich]),
            EstimatedAtNodesScaled=sum(node$mean[,yearIndex]),
            
            ObservedOffset=mean(observedOffset[yearWhich]) / offsetScale,
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
      x <- stat[!colnames(stat) %in% c("Year")]
      print(colSums(x))
      message("Column means:")
      print(colMeans(x))
      message("Correlations:")
      print(cor(x))
      
      return(invisible(stat))
    },
    
    collectHyperparameters = function() {
      library(INLA)
      
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
      if (any(rownames(result$summary.hyperpar)=="GroupRho for st")) # AR1
        x <- rbind(x, rho=result$summary.hyperpar["GroupRho for st",])
      
      return(x)
    },
    
    samplePosterior = function(n, index) {
      indexObserved <- inla.stack.index(fullStack, "observed")$data
      return(callSuper(n=n, index=indexObserved))
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
    
    getPopulationDensityInterpolate.internal = function(xyzt, templateRaster=study$getTemplateRaster(), maskPolygon=study$studyArea$boundary) {
      if (any(is.infinite(xyzt$z)) | any(is.nan(xyzt$z))) return(NULL)
      
      cellArea <- prod(res(templateRaster)) # m^2
      densityRaster <- SpatioTemporalRaster(study=study)$interpolate(xyzt, templateRaster=templateRaster, transform=sqrt, inverseTransform=square, boundary=maskPolygon, layerNames=sort(unique(xyzt[,4])), weights=cellArea)
      return(densityRaster)
    },
    
    getPopulationDensityInterpolate = function(templateRaster=study$getTemplateRaster(), maskPolygon=study$studyArea$boundary, getSD=TRUE) {
      if (is.null(data$fittedMean))
        stop("Did you forgot to run collectEstimates() first?")
      
      library(raster)
      
      xyztMean <- data.frame(getUnscaledObservationCoordinates(), z=data$fittedMean / offsetScale, t=data$year)
      meanPopulationDensityRaster <- getPopulationDensityInterpolate.internal(xyztMean, templateRaster=templateRaster, maskPolygon=maskPolygon)
      sdPopulationDensityRaster <- if (getSD) {
        xyztSD <- data.frame(getUnscaledObservationCoordinates(), z=data$fittedSD / offsetScale, t=data$year)
        getPopulationDensityInterpolate.internal(xyztSD, templateRaster=templateRaster, maskPolygon=maskPolygon)
      }
      else SpatioTemporalRaster$new(study=study)
      
      return(invisible(list(mean=meanPopulationDensityRaster, sd=sdPopulationDensityRaster)))
    },
    
    # Unreliable?!?!
    getPopulationDensityAtMesh = function(templateRaster=study$getTemplateRaster(), maskPolygon=study$studyArea$boundary, getSD=TRUE) {
      if (length(node) == 0)
        stop("Did you forgot to run collectEstimates() first?")
      
      library(raster)

      meanPopulationDensityRaster <- SpatioTemporalRaster$new(study=study)
      sdPopulationDensityRaster <- SpatioTemporalRaster$new(study=study)
      cellArea <- prod(res(templateRaster)) # m^2
      
      for (year in sort(unique(data$year))) {
        yearIndex <- year - min(data$year) + 1
        message("Processing year ", year, "...")
        
        meanRaster <- project(projectValues=node$mean[,yearIndex] / offsetScale, projectionRaster=templateRaster, maskPolygon=maskPolygon, weights=cellArea)
        meanPopulationDensityRaster$addLayer(meanRaster, year)
        
        if (getSD) {
          sdRaster <- project(projectValues=node$sd[,yearIndex] / offsetScale, projectionRaster=templateRaster, maskPolygon=maskPolygon, weights=cellArea)
          sdPopulationDensityRaster$addLayer(sdRaster, year)
        }
      }
      
      return(invisible(list(mean=meanPopulationDensityRaster, sd=sdPopulationDensityRaster)))
    }
  )
)

SimulatedFMPModel <- setRefClass(
  Class = "SimulatedFMPModel",
  contains = "FMPModel",
  fields = list(
    iteration = "integer"
  ),
  methods = list(
    getEstimatesFileIterations = function() {
      if (inherits(study, "undefinedField") | length(modelName) == 0)
        stop("Provide study and modelName parameters.")
      return(study$context$getIterationIds(dir=study$context$scratchDirectory, name=modelName, response=study$response, region=study$studyArea$region, tag="(\\d+)"))
    },
    
    getEstimatesFileName = function() {
      if (inherits(study, "undefinedField") | length(modelName) == 0 | length(iteration) == 0)
        stop("Provide study, modelName and iteration parameters.")
      return(study$context$getLongFileName(study$context$scratchDirectory, name=modelName, response=study$response, region=study$studyArea$region, tag=iteration))
    },
    
    samplePosterior = function(n, index) {
      stop("Unsupported.")
    }
  )
)

SimulatedSmoothModelMeanTemporal <- setRefClass(
  Class = "SimulatedSmoothModelMeanTemporal",
  contains = "SmoothModelMeanTemporal",
  fields = list(
    iteration = "integer"  
  ),
  methods = list(
    getEstimatesFileIterations = function() {
      if (inherits(study, "undefinedField") | length(modelName) == 0)
        stop("Provide study and modelName parameters.")
      return(study$context$getIterationIds(dir=study$context$scratchDirectory, name=modelName, response=study$response, region=study$studyArea$region, tag="(\\d+)"))
    },
    
    getEstimatesFileName = function() {
      if (inherits(study, "undefinedField") | length(modelName) == 0 | length(iteration) == 0)
        stop("Provide study, modelName and iteration parameters.")
      return(study$context$getLongFileName(study$context$scratchDirectory, name=modelName, response=study$response, region=study$studyArea$region, tag=iteration))
    }
  )
)

SimulatedSmoothModelTemporal <- setRefClass(
  Class = "SimulatedSmoothModelTemporal",
  contains = "SmoothModelTemporal",
  fields = list(
    iteration = "integer"
  ),
  methods = list(
    getEstimatesFileIterations = function() {
      if (inherits(study, "undefinedField") | length(modelName) == 0)
        stop("Provide study and modelName parameters.")
      return(study$context$getIterationIds(dir=study$context$scratchDirectory, name=modelName, response=study$response, region=study$studyArea$region, tag="(\\d+)"))
    },
    
    getEstimatesFileName = function() {
      if (inherits(study, "undefinedField") | length(modelName) == 0 | length(iteration) == 0)
        stop("Provide study, modelName and iteration parameters.")
      return(study$context$getLongFileName(study$context$scratchDirectory, name=modelName, response=study$response, region=study$studyArea$region, tag=iteration))
    }
  )
)

SimulatedSmoothModelSpatioTemporal <- setRefClass(
  Class = "SimulatedSmoothModelSpatioTemporal",
  contains = "SmoothModelSpatioTemporal",
  fields = list(
    iteration = "integer"
  ),
  methods = list(
    getEstimatesFileIterations = function() {
      if (inherits(study, "undefinedField") | length(modelName) == 0)
        stop("Provide study and modelName parameters.")
      return(study$context$getIterationIds(dir=study$context$scratchDirectory, name=modelName, response=study$response, region=study$studyArea$region, tag="(\\d+)"))
    },
    
    getEstimatesFileName = function() {
      if (inherits(study, "undefinedField") | length(modelName) == 0 | length(iteration) == 0)
        stop("Provide study, modelName and iteration parameters.")
      return(study$context$getLongFileName(study$context$scratchDirectory, name=modelName, response=study$response, region=study$studyArea$region, tag=iteration))
    }
  )
)

FinlandFMPModel <- setRefClass(
  Class = "FinlandFMPModel",
  contains = c("FMPModel", "FinlandCovariates"),
  fields = list(
  ),
  methods = list(
    initialize = function(...) {
      callSuper(covariatesName="FinlandSmoothModelCovariates", ...)
      return(invisible(.self))
    },
    
    predictDistances = function(formula=study$getDistanceCovariatesModel(), intervalH=study$getTrackSampleInterval()) {
      distances <- study$predictDistances(formula=formula, data=covariates, intervalH=intervalH)
      return(distances)
    }
  )
)

FinlandSmoothModelTemporal <- setRefClass(
  Class = "FinlandSmoothModelTemporal",
  contains = c("SmoothModelTemporal", "FinlandCovariates"),
  fields = list(
  ),
  methods = list(
    initialize = function(...) {
      callSuper(covariatesName="FinlandSmoothModelCovariates", ...)
      return(invisible(.self))
    },
    
    predictDistances = function(formula=study$getDistanceCovariatesModel(), intervalH=study$getTrackSampleInterval()) {
      distances <- study$predictDistances(formula=formula, data=covariates, intervalH=intervalH)
      return(distances)
    }    
  )
)

FinlandSmoothModelSpatioTemporal <- setRefClass(
  Class = "FinlandSmoothModelSpatioTemporal",
  contains = c("SmoothModelSpatioTemporal", "FinlandCovariates"),
  fields = list(
  ),
  methods = list(
    initialize = function(...) {
      callSuper(covariatesName="FinlandSmoothModelCovariates", ...)
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

FinlandSmoothModelSpatioTemporalNoDistances <- setRefClass(
  Class = "FinlandSmoothModelSpatioTemporalNoDistances",
  contains = c("SmoothModelSpatioTemporal", "FinlandCovariates"),
  fields = list(
  ),
  methods = list(
    initialize = function(...) {
      callSuper(covariatesName="FinlandSmoothModelCovariates", ...)
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

FinlandRussiaSmoothModelSpatioTemporal <- setRefClass(
  Class = "FinlandRussiaSmoothModelSpatioTemporal",
  contains = c("SmoothModelSpatioTemporal"),
  fields = list(
  ),
  methods = list(
  )
)

##################################################################################
# OLD STUFF
if (F) {
  
  
  TODO_getPopulationSize = function(populationDensity, tracks, habitatWeights) {
    if (is.null(data$fittedMean))
      stop("Did you forgot to run collectEstimates() first?")
    
    if (!missing(habitatWeights)) if (!is.null(habitatWeights)) populationDensity$weight(habitatWeights)
    
    # TODO: this can be done faster without the rasters, but need the area of the study area
    #area <- maskPolygon@polygons[[1]]@area
    #if (!missing(habitatWeights)) if (!is.null(habitatWeights)) area <- cellStats(habitatWeights, sum)      
    #meanPopulationSize <- ddply(data, .(year), function(x, area)
    #  data.frame(z=mean(x$fittedMean / offsetScale) * area, year=x$year[1], area))
    #populationSize <- SimulationPopulationSize(study=study, iteration=iteration, modelName=modelName)
    
    populationSize <- populationDensity$integrate(volume=SimulationPopulationSize$new(study=study, iteration=iteration, modelName=modelName))
    if (missing(tracks)) populationSize$loadValidationData()
    else populationSize$loadValidationData(tracks=tracks)
    
    return(invisible(populationSize))
  }
  
  
  #getPopulationSize = function(tracks, withHabitatWeights=FALSE) {
  #  if (length(node) == 0)
  #    stop("Run collectEstimates() first.")
  #  populationDensity <- getPopulationDensity(getSD=FALSE)
  #  
  #  if (withHabitatWeights) {
  #    habitatWeights <- CORINEHabitatWeights$new(study=study)
  #    if (missing(tracks)) tracks <- study$loadTracks(iteration=iteration)
  #    habitatSelection <- tracks$getHabitatPreferences(habitatWeightsTemplate=habitatWeights, nSamples=30, save=FALSE) # TODO: save
  #    habitatWeights$setHabitatSelectionWeights(habitatSelection)
  #    habitatWeightsRaster <- habitatWeights$getWeightsRaster(save=FALSE)
  #    populationDensity$mean$weight(habitatWeightsRaster)
  #  }
  #  
  #  populationSize <- populationDensity$mean$integrate(volume=SimulationPopulationSize$new(study=study, iteration=iteration, modelName=modelName))
  #  if (missing(tracks)) populationSize$loadValidationData()
  #  else populationSize$loadValidationData(tracks=tracks)
  # 
  #  return(invisible(populationSize))
  #}
  
  
  
  #if (withHabitatWeights) {
  #  habitatWeights <- CORINEHabitatWeights$new(study=study)
  #  if (missing(tracks)) tracks <- study$loadTracks(iteration=iteration)
  #  habitatSelection <- tracks$getHabitatPreferences(habitatWeightsTemplate=habitatWeights, nSamples=30, save=FALSE) # TODO: save
  #  habitatWeights$setHabitatSelectionWeights(habitatSelection)
  #  habitatWeightsRaster <- habitatWeights$getWeightsRaster(save=FALSE)
  #  populationDensity$mean$weight(habitatWeightsRaster)
  #}
  
  
  switchToMesh = function() {
    assign("getPopulationDensityAtObserved", .self$getPopulationDensity, envir=as.environment(.self))
    assign("getPopulationDensity", .self$getPopulationDensityAtMesh, envir=as.environment(.self))
  }
  
  switchToObservations = function() {
    assign("getPopulationDensityAtMesh", .self$getPopulationDensity, envir=as.environment(.self))
    assign("getPopulationDensity", .self$getPopulationDensityAtObserved, envir=as.environment(.self))
  }
  
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
  }
  
  saveMeshNodeCovariates = function(save=FALSE) {
    stop("Method saveMeshNodeCovariates unimplemented.")
  }
  
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
  }
  
  plotProjection = function(projectionRaster=study$getTemplateRaster(), variable, weights=1) {
    meanRaster <- project(projectValues=variable, projectionRaster=projectionRaster, weights=weights)
    op <- par(mar=rep(0, 4))
    plot(meanRaster)
    plot(study$studyArea$boundary, add=T)
    points(locations / coordsScale, col=data$intersections)
    par(op)
  }
  
  plotSpatialEstimatedMean = function(yearIndex=1) {
    plotProjection(variable=node$mean[,yearIndex], weights=prod(res(study$getTemplateRaster())))
  }
  
  plotSpatialEstimatedSD = function(yearIndex=1) {
    plotProjection(variable=node$sd[,yearIndex], weights=prod(res(study$getTemplateRaster()))^2)
  }
  
  plotSpatialRandomEffectMean = function(yearIndex=1) {
    plotProjection(variable=node$spatialMean[,yearIndex])
  }
  
  plotSpatialRandomEffectSD = function(yearIndex=1) {
    plotProjection(variable=node$spatialSd[,yearIndex])
  }
  
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
}
