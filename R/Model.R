library(INLA)
setOldClass("inla.mesh")
setOldClass("inla.spde2")
setOldClass("inla.data.stack")
setOldClass("inla")

Model <- setRefClass(
  Class = "Model",
  fields = list(
    name = "character",
    study = "Study",
    data = "data.frame",
    locations = "matrix",
    covariates = "data.frame",
    model = "formula",
    coordsScale = "numeric",
    years = "integer",
    mesh = "inla.mesh",
    spde = "inla.spde2",
    index = "list",
    A = "Matrix",
    dataStack = "inla.data.stack",
    offset = "numeric",
    result = "inla"
  ),
  methods = list(
    initialize = function(...) {
      callSuper(...)
      return(invisible(.self))
    },
    
    getResultFileName = function() {
      if (inherits(study, "undefinedField") | length(name) == 0)
        stop("Provide study and name parameters.")
      return(study$context$getLongFileName(dir=study$context$resultDataDirectory, name="Result", response=response, region=study$studyArea$region, tag=name))
    },
    
    saveResult = function() {
      save(name, locations, covariates, model, coordsScale, years, mesh, spde, index, A, dataStack, offset, result, file=getResultFileName())
    },
    
    loadResult = function() {
      load(getResultFileName(), envir=as.environment(.self))
    }
  )
)

SmoothModel <- setRefClass(
  Class = "SmoothModel",
  contains = "Model",
  fields = list(
  ),
  methods = list(
    initialize = function(...) {
      callSuper(name="Smooth", ...)
      return(invisible(.self))
    },
    
    setup = function(intersections, meshParams) {
      library(INLA)
            
      message("Constructing mesh...")      

      data <<- subset(intersections$intersections, response == study$response)@data
      coordsScale <<- meshParams$coordsScale
      locations <<- coordinates(intersections$intersections) * coordsScale
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
      dataStack <<- inla.stack(data=list(response=as.numeric(data$response)),
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
      plot(surveyRoutes$surveyRoutes, col="blue", add=T)
    },

    estimate = function(family="nbinomial") {
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
    },
    
    collectResults = function() {
      library(INLA)
      library(plyr)
      
      indexObserved <- inla.stack.index(dataStack, "observed")$data
      #fitted <- exp(result$summary.linear.predictor$mean[indexObserved]) # exp(E(eta))
      eta <- result$marginals.linear.predictor[indexObserved]
      fitted <- ldply(eta, function(eta) getINLAEstimates(eta, exp))[,-1] # E(exp(eta))
      data$eta <<- result$summary.linear.predictor$mean[indexObserved]
      data$fittedMean <<- fitted$mean
      data$fittedSD <<- fitted$sd
      
      spdeResult <- inla.spde2.result(result, "st", spde)
      logKappa <- getINLAEstimates(spdeResult$marginals.log.kappa[[1]], coordsScale=1)
      logTau <- getINLAEstimates(spdeResult$marginals.log.tau[[1]], coordsScale=1)
      sd <- getINLAEstimates(spdeResult$marginals.log.variance.nominal[[1]], fun=function(x) sqrt(exp(x)))
      range <- getINLAEstimates(spdeResult$marginals.range.nominal[[1]], coordsScale=coordsScale)
      # sqrt(8)/kappa
      # 1/(4*pi*kappa^2*tau^2)
      
      x <- rbind(kappa=logKappa,
                 tau=logTau,
                 sd=sd,
                 range=range)
      if (any(rownames(result$summary.hyperpar)=="GroupRho for st"))
        x <- rbind(x, rho=result$summary.hyperpar["GroupRho for st",])

      return(invisible(x))
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
      
      if (!missing(maskPolygon)) projectionRaster <- mask(projectionRaster, maskPolygon)
      
      return(invisible(projectionRaster))
    },
    
    getPopulationDensity = function(projectionRaster) {
      library(raster)

      populationDensityRaster <- SpatioTemporalRaster$new()
      
      intersectionsMean <- matrix(data=data$fittedMean, nrow=length(years), ncol=mesh$n)
      intersectionsSD <- matrix(data=data$fittedSD, nrow=length(years), ncol=mesh$n)
      
      for (yearIndex in 1:nrow(intersectionsMean)) {
        year <- as.character(yearIndex + years[1] - 1)
        message("Processing year ", year, "...")
        
        meanRasterLayer <- project(estimates=intersectionsMean[yearIndex,], projectionRaster=projectionRaster, maskPolygon=study$studyArea$boundary)
        names(meanRasterLayer) <- year
        #populationDensityMean[[yearIndex]] <- x
        
        sdRasterLayer <- project(estimates=intersectionsSD[yearIndex,], projectionRaster=projectionRaster, maskPolygon=study$studyArea$boundary)
        names(sdRasterLayer) <- year
        #populationDensitySD[[yearIndex]] <- x
        
        populationDensityRaster$addRasters(meanRasterLayer=meanRasterLayer, sdRasterLayer=sdRasterLayer)
      }
      
      #populationDensityRaster$setRasters(populationDensityMean, populationDensitySD)
      return(invisible(populationDensityRaster))
    }    
  )
)

FinlandSmoothModel <- setRefClass(
  Class = "FinlandSmoothModel",
  contains = "SmoothModel",
  fields = list(
  ),
  methods = list(
  )
)
