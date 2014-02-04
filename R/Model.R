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
    locations = "matrix",
    covariates = "data.frame",
    model = "formula",
    coordsScale = "numeric",
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
      return(study$context$getLongFileName(dir=study$context$resultDataDirectory, name="Result", response=response, tag=name, region=study$studyArea$region))
    },
    
    saveResult = function() {
      save(name, locations, covariates, model, coordsScale, mesh, spde, index, A, dataStack, offset, result, file=getResultFileName())
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

      data <- subset(intersections$intersections, response == study$response)
            
      message("Constructing mesh...")      
      
      coordsScale <<- meshParams$coordsScale
      locations <<- coordinates(data) * coordsScale
      mesh <<- inla.mesh.create.helper(points.domain=locations,
                                       min.angle=meshParams$minAngle,
                                       max.edge=meshParams$maxEdge * coordsScale,
                                       cutoff=meshParams$cutOff * coordsScale)
      
      message("Number of nodes in the mesh = ", mesh$n)
      
      spde <<- inla.spde2.matern(mesh)
      model <<- response ~ -1 + intercept + f(st, model=spde, group=st.group, control.group=list(model="ar1"))

      years <- data$year
      nYears <- length(unique(years))
      groupYears <- as.integer(years - min(years) + 1)
      
      index <<- inla.spde.make.index("st", n.spde=mesh$n, n.group=nYears)
      A <<- inla.spde.make.A(mesh, loc=locations, group=groupYears, n.group=nYears)
      dataStack <<- inla.stack(data=list(response=as.numeric(data$response)),
                               A=list(A),
                               effects=list(c(index, list(intercept=1))),
                               tag="observed")
      
      offset <<- 2/pi * data$length * data$duration * data$distance
      
      return(invisible(.self))
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
        stop("INLA failed to run.")
      }

      library(plyr)
      indexObserved <- inla.stack.index(dataStack, "observed")$data
      #fitted <- exp(result$summary.linear.predictor$mean[indexObserved]) # TODO: We want E(exp(eta)), not exp(E(eta)) !!!
      eta <- result$marginals.linear.predictor[indexObserved]
      fitted <- laply(eta, function(eta) inla.emarginal(function(x) exp(x), eta))
      return(fitted)
    },
    
    plotMesh = function(surveyRoutes) {
      library(INLA)
      library(sp)
      meshUnscaled <- mesh
      meshUnscaled$loc <- mesh$loc / coordsScale
      plot(meshUnscaled, asp=1, main="")
      plot(study$studyArea$boundary, border="darkgray", add=T)
      plot(surveyRoutes$surveyRoutes, col="blue", add=T)
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
