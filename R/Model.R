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
    dataStack = "inla.data.stack",
    offset = "numeric",
    result = "inla",
    node = "list"
  ),
  methods = list(
    initialize = function(...) {
      callSuper(...)
      return(invisible(.self))
    },
    
    #getResultFileName = function() {
    #  if (inherits(study, "undefinedField") | length(name) == 0)
    #    stop("Provide study and name parameters.")
    #  return(study$context$getLongFileName(dir=study$context$resultDataDirectory, name="Result", response=study$response, region=study$studyArea$region, tag=name))
    #},
    
    saveResult = function(fileName) {
      #fileName <- getItemFileName(directory, name)
      #fileName <- getResultFileName()
      message("Saving result to ", fileName, "...")
      save(locations, data, model, coordsScale, years, mesh, spde, index, A, dataStack, offset, result, file=fileName)
    },
    
    loadResult = function(fileName) {
      load(fileName, envir=as.environment(.self))
      #load(getResultFileName(), envir=as.environment(.self))
      return(invisible(.self))
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
      callSuper(...)
      return(invisible(.self))
    },
    
    setup = function(intersections, meshParams, family="nbinomial") {
      library(INLA)
      
      family <<- family
      
      message("Constructing mesh...")

      intersectionsSubset <- subset(intersections$intersections, response == study$response)
      data <<- intersectionsSubset@data
      coordsScale <<- meshParams$coordsScale
      locations <<- coordinates(intersectionsSubset) * coordsScale
      mesh <<- inla.mesh.create.helper(points.domain=locations,
                                       min.angle=meshParams$minAngle,
                                       max.edge=meshParams$maxEdge * coordsScale,
                                       cutoff=meshParams$cutOff * coordsScale)
      
      message("Number of nodes in the mesh = ", mesh$n)
      
      spde <<- inla.spde2.matern(mesh)
      model <<- response ~ -1 + intercept + f(st, model=spde, group=st.group, control.group=list(model="ar1"))

      years <<- as.integer(sort(unique(data$year)))
      nYears <- length(years)
      groupYears <- as.integer(data$year - min(years) + 1)
      
      index <<- inla.spde.make.index("st", n.spde=mesh$n, n.group=nYears)
      A <<- inla.spde.make.A(mesh, loc=locations, group=groupYears, n.group=nYears)
      dataStack <<- inla.stack(data=list(response=data$intersections),
                               A=list(A),
                               effects=list(c(index, list(intercept=1))),
                               tag="observed")
      
      offset <<- 2/pi * data$length * data$duration * data$distance
      
      return(invisible(.self))
    },
    
    plotMesh = function(surveyRoutes) {
      library(INLA)
      library(sp)
      meshUnscaled <- mesh
      meshUnscaled$loc <- mesh$loc / coordsScale
      plot(meshUnscaled, asp=1, main="")
      plot(study$studyArea$boundary, border="darkgray", add=T)
      
      if (missing(surveyRoutes)) {}
        # TODO: plot intersection coordinates
      else plot(surveyRoutes$surveyRoutes, col="blue", add=T)
      
      return(invisible(.self))
    },

    estimate = function(save=FALSE, fileName) {
      library(INLA)
      
      result <<- inla(model,
                      family=family,
                      data=inla.stack.data(dataStack),
                      E=.self$offset,
                      verbose=TRUE,
                      control.predictor=list(A=inla.stack.A(dataStack), compute=TRUE),
                      control.compute=list(cpo=FALSE, dic=TRUE))
      
      if (is.null(result$ok) || result$ok == FALSE) {
        warning("INLA failed to run.")
      }
      else {
        if (save) saveResult(fileName=fileName)
      }
    },

    # TODO: fix this if needed
    collectResults = function(quick=FALSE) {
      if (quick) return(collectResultsQuick())
      
      library(INLA)
      library(plyr)
      
      message("Processing fitted values...")
      indexObserved <- inla.stack.index(dataStack, "observed")$data

      data$eta <<- result$summary.linear.predictor$mean[indexObserved]
      data$fittedMean <<- result$summary.fitted.values$mean[indexObserved]
      data$fittedSD <<- result$summary.fitted.values$sd[indexObserved]
      
      #message("Processing random effects...")
      ## TODO: this is wrong, fix!
      #intercept <- result$summary.fixed["intercept","mean"]
      #st <- ldply(result$marginals.random$st,
      #            function(st) getINLAEstimates(st, function(st) exp(st + intercept)), .parallel=TRUE)[,-1]
      #node <<- list()
      #node$mean <<- matrix(data=st$mean, nrow=length(years), ncol=mesh$n)
      #node$sd <<- matrix(data=st$sd, nrow=length(years), ncol=mesh$n)
      
      
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
    },

    collectResultsQuick = function() {
      library(INLA)
      
      indexObserved <- inla.stack.index(dataStack, "observed")$data
      data$eta <<- result$summary.linear.predictor$mean[indexObserved]
      data$fittedMean <<- result$summary.fitted.values$mean[indexObserved]
      data$fittedSD <<- result$summary.fitted.values$sd[indexObserved]
      
      node <<- list()
      node$mean <<- matrix(data=st$mean, nrow=length(years), ncol=mesh$n)
      node$sd <<- matrix(data=st$sd, nrow=length(years), ncol=mesh$n)
      
      spdeResult <- inla.spde2.result(result, "st", spde)
      logKappa <- spdeResult$summary.log.kappa
      logTau <- spdeResult$summary.log.tau
      sd <- sqrt(exp(spdeResult$summary.log.variance.nominal))
      range <- exp(spdeResult$summary.log.range.nominal) * coordsScale
      x <- rbind(logKappa=logKappa,
                 logTau=logTau,
                 sd=sd,
                 range=range)
      if (any(rownames(result$summary.hyperpar)=="GroupRho for st"))
        x <- rbind(x, rho=result$summary.hyperpar["GroupRho for st",])
      return(x)      
    },
    
    project = function(estimates, projectionRaster, maskPolygon) {
      library(INLA)
      library(raster)
      
      projector <- inla.mesh.projector(mesh,
                                       dims=c(ncol(projectionRaster), nrow(projectionRaster)),
                                       xlim=c(xmin(projectionRaster), xmax(projectionRaster)),
                                       ylim=c(ymin(projectionRaster), ymax(projectionRaster)))
      projectedEstimates <- inla.mesh.project(projector, estimates)
      area <- prod(res(projectionRaster))
      values(projectionRaster) <- t(projectedEstimates[,ncol(projectedEstimates):1]) * area
      
      if (!missing(maskPolygon) & !is.null(maskPolygon))
        projectionRaster <- mask(projectionRaster, maskPolygon)
      
      return(invisible(projectionRaster))
    },
    
    getPopulationDensity = function(projectionRaster, maskPolygon=study$studyArea$boundary, getSD=TRUE) {
      if (length(node) == 0)
        stop("Did you forgot to run collectResults() first?")
      
      library(raster)

      meanPopulationDensityRaster <- SpatioTemporalRaster$new()
      sdPopulationDensityRaster <- SpatioTemporalRaster$new()      
      
      for (yearIndex in 1:nrow(node$mean)) {
        year <- as.character(yearIndex + years[1] - 1)
        message("Processing year ", year, "...")
        
        meanRaster <- project(estimates=node$mean[yearIndex,], projectionRaster=projectionRaster, maskPolygon=maskPolygon)
        names(meanRaster) <- year
        meanPopulationDensityRaster$addLayer(meanRaster)
        
        if (getSD) {
          sdRaster <- project(estimates=node$sd[yearIndex,], projectionRaster=projectionRaster, maskPolygon=maskPolygon)
          names(sdRaster) <- year
          sdPopulationDensityRaster$addLayer(sdRaster)
        }
      }
      
      return(invisible(list(mean=meanPopulationDensityRaster, sd=sdPopulationDensityRaster)))
    }    
  )
)
