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
    modelName = "character"
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
      save(locations, data, model, coordsScale, years, mesh, spde, index, A, fullStack, result, file=fileName)
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
    }
  )
)

SmoothModel <- setRefClass(
  Class = "SmoothModel",
  contains = "Model",
  fields = list(
    family = "character"
  ),
  methods = list(
    initialize = function(...) {
      callSuper(modelName="SmoothModel", ...)
      return(invisible(.self))
    },
    
    setup = function(intersections, meshParams, family="nbinomial") {
      library(INLA)
      
      message("Constructing mesh...")

      data <<- intersections$getData()
      coordsScale <<- meshParams$coordsScale
      locations <<- intersections$getCoordinates() * coordsScale
      mesh <<- inla.mesh.create.helper(points.domain=locations,
                                       min.angle=meshParams$minAngle,
                                       max.edge=meshParams$maxEdge * coordsScale,
                                       cutoff=meshParams$cutOff * coordsScale)
      
      message("Number of nodes in the mesh = ", mesh$n)
      
      family <<- family
      model <<- response ~ -1 + intercept + f(st, model=spde, group=st.group, control.group=list(model="ar1"))

      years <<- as.integer(sort(unique(data$year)))
      nYears <- length(years)
      groupYears <- as.integer(data$year - min(years) + 1)
      
      spde <<- inla.spde2.matern(mesh)
      index <<- inla.spde.make.index("st", n.spde=mesh$n, n.group=nYears)
      
      A <<- inla.spde.make.A(mesh, loc=locations, group=groupYears, n.group=nYears)
      obsStack <<- inla.stack(data=list(response=data$intersections,
                                        E=2/pi * data$length * data$duration * data$distance,
                                        link=1),
                               A=list(A),
                               effects=list(c(index, list(intercept=1))),
                               tag="observed")
      
      predStack <<- inla.stack(data=list(response=NA,
                                         E=rep(2/pi * 12000 * 1 * mean(data$distance), mesh$n * nYears),
                                         link=1),
                               A=list(1),
                               effects=list(c(index, list(intercept=1))),
                               tag="predicted")
            
      return(invisible(.self))
    },
    
    estimate = function(save=FALSE, fileName=getEstimatesFileName()) {
      library(INLA)
      
      fullStack <<- inla.stack(obsStack, predStack)
      stackData <- inla.stack.data(fullStack)
      result <<- inla(model,
                      family=family,
                      data=inla.stack.data(fullStack, spde=spde),
                      E=stackData$E,
                      verbose=TRUE,
                      control.predictor=list(A=inla.stack.A(fullStack), link=stackData$link, compute=TRUE),
                      control.compute=list(cpo=FALSE, dic=TRUE))
      
      if (is.null(result$ok) || result$ok == FALSE) {
        warning("INLA failed to run.")
      }
      else {
        if (save) saveEstimates(fileName=fileName)
      }
    },

    collectEstimates = function(weightsAtSurveyRoutes=1, weightsAtNodes=1, quick=FALSE) {
      library(INLA)
      library(plyr)
      
      message("Processing fitted values...")
      
      indexObserved <- inla.stack.index(fullStack, "observed")$data
      data$eta <<- result$summary.linear.predictor$mean[indexObserved] - log(weightsAtSurveyRoutes)
      data$fittedMean <<- result$summary.fitted.values$mean[indexObserved] / weightsAtSurveyRoutes
      data$fittedSD <<- result$summary.fitted.values$sd[indexObserved] / weightsAtSurveyRoutes^2
      
      message("Processing random effects...")
        
      yearsVector <- sort(unique(data$year))
      indexPredicted <- inla.stack.index(fullStack, "predicted")$data
      node$mean <<- matrix(result$summary.fitted.values$mean[indexPredicted] / weightsAtNodes, nrow=mesh$n, ncol=length(yearsVector))
      node$sd <<- matrix(result$summary.fitted.values$sd[indexPredicted] / weightsAtNodes^2, nrow=mesh$n, ncol=length(yearsVector))
      
      
      a <- inla.emarginal(exp, result$marginals.fixed$intercept)
      if (!all(a >= exp(result$summary.fixed["intercept","mean"])))
        warning("Jensen's inequality does not hold for fixed effects.")
      
      b <- numeric(length(result$marginals.random$st))
      for (i in 1:length(result$marginals.random$st)) {
        b[i] <- inla.emarginal(exp, result$marginals.random$st[[i]])
      }
      if (!all(b >= exp(result$summary.random$st$mean)))
        warning("Jensen's inequality does not hold for random effects.")
      
      e <- a * b
      e.m <- matrix(e / weightsAtNodes, nrow=mesh$n, ncol=length(yearsVector))
      if (!all(node$mean >= e.m))
        warning("Jensen's inequality does not hold for fitted values.")
      
      
      # TODO: fix this
if (F) {      
      message("Processing hyperparameters...")
      spdeResult <- inla.spde2.result(result, "st", spde)
      logKappa <- getINLAEstimates(spdeResult$marginals.log.kappa[[1]], coordsScale=1)
      logTau <- getINLAEstimates(spdeResult$marginals.log.tau[[1]], coordsScale=1) ## scale ??
      sd <- getINLAEstimates(spdeResult$marginals.log.variance.nominal[[1]], fun=function(x) sqrt(exp(x)))
      range <- getINLAEstimates(spdeResult$marginals.range.nominal[[1]], coordsScale=coordsScale)
      # sqrt(8)/kappa
      # 1/(4*pi*kappa^2*tau^2)
      
      x <- rbind(logKappa=logKappa,
                 logTau=logTau,
                 sd=sd,
                 range=range)
      if (any(rownames(result$summary.hyperpar)=="GroupRho for st"))
        x <- rbind(x, rho=result$summary.hyperpar["GroupRho for st",])
            
      return(x)
}
    },

    getPredictedIntersections = function() {
      x <- data
      indexObserved <- inla.stack.index(fullStack, "observed")$data
      x$offset <- inla.stack.data(fullStack)$E[indexObserved]
      intersections <- ddply(x, .(year), function(x) {
        data.frame(Observed=sum(x$intersections), Predicted=sum(x$fittedMean * x$offset), logObsOff=log(sum(x$intersections / x$offset)), eta=sum(x$eta))
      })
      return(intersections)
    },
    
    project = function(projectValues, projectionRaster, maskPolygon) {
      library(INLA)
      library(raster)
      
      projector <- inla.mesh.projector(getUnscaledMesh(),
                                       dims=c(ncol(projectionRaster), nrow(projectionRaster)),
                                       xlim=c(xmin(projectionRaster), xmax(projectionRaster)),
                                       ylim=c(ymin(projectionRaster), ymax(projectionRaster)))
      projectedEstimates <- inla.mesh.project(projector, projectValues)
      
      cellArea <- prod(res(projectionRaster))
      values(projectionRaster) <- t(projectedEstimates[,ncol(projectedEstimates):1]) * cellArea
      
      if (!missing(maskPolygon) & !is.null(maskPolygon))
        projectionRaster <- mask(projectionRaster, maskPolygon)
      
      return(invisible(projectionRaster))
    },
  
    getPopulationDensity = function(templateRaster=study$getTemplateRaster(), maskPolygon=study$studyArea$boundary, getSD=TRUE) {
      if (length(node) == 0)
        stop("Did you forgot to run collectEstimates() first?")
      
      library(raster)
      #library(ST)

      meanPopulationDensityRaster <- SpatioTemporalRaster$new(study=study)
      sdPopulationDensityRaster <- SpatioTemporalRaster$new(study=study)
      
      xyzMean <- data.frame(Year=data$year, locations / coordsScale, z=data$fittedMean)
      xyzSD <- data.frame(Year=data$year, locations / coordsScale, z=data$fittedSD)
      #cellArea <- prod(res(templateRaster)) # m^2
      
      for (year in sort(unique(xyzMean$Year))) {
        message("Processing year ", year, "...")

        yearIndex <- year - min(xyzMean$Year) + 1
        meanRaster <- project(projectValues=node$mean[,yearIndex], projectionRaster=templateRaster, maskPolygon=maskPolygon)
        #meanRaster <- rasterInterpolate(subset(xyzMean, Year == year)[,-1], templateRaster=templateRaster, transform=sqrt, inverseTransform=square) * cellArea        
        #if (!(missing(maskPolygon) | is.null(maskPolygon)))
          #meanRaster <- mask(meanRaster, maskPolygon)
        meanPopulationDensityRaster$addLayer(meanRaster, year)
        
        if (getSD) {
          sdRaster <- project(projectValues=node$sd[,yearIndex], projectionRaster=templateRaster, maskPolygon=maskPolygon)
          #sdRaster <- rasterInterpolate(subset(xyzSD, Year == year)[,-1], templateRaster=templateRaster, transform=sqrt, inverseTransform=square) * cellArea
          #if (!(missing(maskPolygon) | is.null(maskPolygon)))
            #sdRaster <- mask(sdRaster, maskPolygon)
          sdPopulationDensityRaster$addLayer(sdRaster, year)
        }
      }
      
      return(invisible(list(mean=meanPopulationDensityRaster, sd=sdPopulationDensityRaster)))
    },
    
    associateMeshLocationsWithDate = function() {
      library(plyr)
      library(fields)
      
      nYears <- length(years)
      meshLocations <- getUnscaledMesh()$loc[,1:2]
      nMeshLocations <- nrow(meshLocations)
      predictLocations <- repeatMatrix(meshLocations, nYears)      
      id <- rep(1:nMeshLocations, nYears)
      year <- rep(years, each=nMeshLocations)
      
      xyt <- ddply(data, .(year), function(x) {
        message("Finding closest observations for mesh nodes for year ", x$year[1], "...")
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
    
    getPopulationSize = function(withHabitatWeights) {
      collectEstimates()
      populationDensity <- getPopulationDensity(getSD=FALSE)
      
      if (mss$hasHabitatWeights()) {
        habitatWeights <- CORINEHabitatWeights$new(study=study)
        tracks <- study$loadTracks(iteration=iteration)
        habitatSelection <- tracks$getHabitatPreferences(habitatWeightsTemplate=habitatWeights, nSamples=30, save=FALSE) # TODO: save
        habitatWeights$setHabitatSelectionWeights(habitatSelection)
        habitatWeighstRaster <- habitatWeights$getWeightsRaster(save=FALSE)
        populationDensity$mean$weight(habitatWeightsRaster)
      }
      
      populationSize <- populationDensity$mean$integrate(volume=SimulationPopulationSize$new(study=study, iteration=iteration))
      populationSize$loadValidationData()
      
      return(invisible(populationSize))
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
    }
  )
)
