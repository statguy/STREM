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
    
    setup = function(intersections, params, covariatesModel=~1, tag, timeout=30) {
      library(INLA)
      library(R.utils)
      
      callSuper(intersections=intersections, params=params, covariatesModel=covariatesModel, tag=tag)
      setModelName(family=params$family, randomEffect=paste("matern", params$timeModel, sep="-"), tag=tag)
      
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
                               timeout=timeout, onTimeout="error")
      message("Mesh nodes = ", mesh$n)
      
      nYears <- length(years)
      groupYears <- as.integer(data$year - min(years) + 1)
      
      spde <<- if (length(spdeParams) == 0) inla.spde2.matern(mesh)
      else inla.spde2.matern(mesh, B.tau=spdeParams$B.tau, B.kappa=spdeParams$B.kappa, theta.prior.mean=spdeParams$theta.prior.mean, theta.prior.prec=spdeParams$theta.prior.prec, constr=spdeParams$constr)
      index <<- inla.spde.make.index("st", n.spde=spde$n.spde, n.group=nYears)
      

      modelMatrix <- getINLAModelMatrix(covariatesModel, intersections$getData())
      covTerms <- terms(covariatesModel, data=intersections$getData())
      covariates <- colnames(modelMatrix)
      intercept <- if (attr(covTerms, "intercept")[1] == 0) NULL else "intercept"
      randomEffect <- if (length(rhoPrior) == 0) "f(st, model=spde, group=st.group, control.group=list(model='ar1'))"
      else "f(st, model=spde, group=st.group, control.group=list(model='ar1', hyper=rhoPrior))"
      model <<- reformulate(termlabels=c(intercept, covariates, randomEffect), response="response", intercept=FALSE)
      #model <<- if (length(rhoPrior) == 0) response ~ -1 + intercept + f(st, model=spde, group=st.group, control.group=list(model="ar1"))
      #else response ~ -1 + intercept + f(st, model=spde, group=st.group, control.group=list(model="ar1", hyper=rhoPrior))
      A <<- inla.spde.make.A(mesh, loc=locations, group=groupYears, n.group=nYears)
      
      message("Number of nodes in the mesh = ", mesh$n)
      message("Average survey route length = ", mean(data$length))
      message("Average count duration = ", mean(data$duration))
      message("Average animal track length = ", mean(data$distance))
      
      effects <- if (!is.null(intercept)) list(c(index, list(intercept=1))) else list(index)
      AList <- if (!is.null(modelMatrix)) {
        effects[[2]] <- modelMatrix
        list(A, 1)
      }
      else list(A)
      obsStack <<- inla.stack(data=list(response=data$intersections,
                                        E=getObservedOffset(),
                                        link=1),
                              A=AList,
                              effects=effects,
                              tag="observed")
      
      #boundarySamples <- study$studyArea$sampleBoundary()
      #if (!is.null(boundarySamples)) {
      #  coords <- repeatMatrix(coordinates(boundarySamples), length(years)) * coordsScale
      #  yearsData <- rep(years, each=length(boundarySamples))      
      #  groupYears <- as.integer(yearsData - min(yearsData) + 1)
      #  augA <<- inla.spde.make.A(mesh, loc=coords, group=groupYears, n.group=nYears)
      #  augStack <<- inla.stack(data=list(response=NA,
      #                                    E=1,
      #                                    link=1),
      #                          A=list(augA),
      #                          effects=list(c(index, list(intercept=1))),
      #                          tag="augmented")
      #}
      
      #predStack <<- inla.stack(data=list(response=NA,
      #                                   E=1,
      #                                   link=1),
      #                         A=list(1),
      #                         effects=list(c(index, list(intercept=1))),
      #                         tag="predicted")
      
      return(invisible(.self))
    },
    
    estimate = function(save=FALSE, fileName=getEstimatesFileName(), verbose=TRUE) {
      library(INLA)
      
      #fullStack <<- if (is.null(augStack)) inla.stack(obsStack, predStack)
      #else inla.stack(obsStack, augStack, predStack)
      fullStack <<- inla.stack(obsStack)
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
                      control.inla=list(int.strategy="grid"),
                      control.compute=list(waic=TRUE, config=TRUE))
      
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
      
      return(invisible(.self))
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
    
    getPopulationDensity = function(templateRaster=study$getTemplateRaster(), maskPolygon=study$studyArea$boundary, habitatWeights=NULL, index, year, .parallel=TRUE) {
      if (is.null(data$fittedMean))
        stop("Did you forgot to run collectEstimates() first?")
      library(raster)
      
      effortWeights <- if (!is.null(habitatWeights)) {
        if (inherits(study$surveyRoutes, "uninitializedField"))
          stop("You specified habitat weights but survey routes are not available for effort weighting.")
        weightedLengths <- study$surveyRoutes$getWeightedLengths(habitatWeights)
        getLengthWeights(weightedLengths, study$surveyRoutes$lengths)
      }
      else 1
      
      populationDensity <- getDensityEstimates(weights=1/effortWeights, aggregate=FALSE, index=index, year=year)
      cellArea <- prod(res(templateRaster))
      populationDensityRaster <- SpatioTemporalRaster(study=study)$interpolate(populationDensity, templateRaster=templateRaster, transform=sqrt, inverseTransform=square, layerNames=sort(unique(populationDensity$year)), weights=cellArea, .parallel=.parallel)
      #populationDensityRaster <- SpatioTemporalRaster(study=study)$interpolate(populationDensity, templateRaster=templateRaster, transform=sqrt, inverseTransform=square, boundary=maskPolygon, layerNames=sort(unique(populationDensity$year)), weights=cellArea, .parallel=.parallel)
      
      return(invisible(populationDensityRaster))
    },
    
    
    #######
    
    
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
    },
    
    getPopulationSize = function(populationDensity, habitatWeightsRaster=NULL) {
      if (missing(populationDensity))
        stop("Required argument 'populationDensity' missing.")
      if (!inherits(populationDensity, "SpatioTemporalRaster"))
        stop("Argument 'populationDensity' must be of type 'SpatioTemporalRaster'")
      if (!is.null(habitatWeightsRaster)) populationDensity$weight(habitatWeightsRaster)
      
      populationSize <- populationDensity$integrate(volume=SimulationPopulationSize(study=study, iteration=iteration, modelName=modelName))
      populationSize$loadValidationData()
      
      return(invisible(populationSize))
    }
  )
)
