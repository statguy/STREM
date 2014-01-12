library(INLA)
setOldClass("inla.mesh")
setOldClass("inla.spde2")

InitialPopulation <- setRefClass(
  Class = "InitialPopulation",
  fields = list(
  ),
  methods = list(
    randomize = function(n) {
      stop("randomize() undefined.")
    },
    
    plotLocations = function(locations) {
      stop("plotLocations() undefined.")
    }
  )
)

RandomInitialPopulation <- setRefClass(
  Class = "RandomInitialPopulation",
  contains = "InitialPopulation",
  fields = list(
    boundary = "SpatialPolygons"
  ),
  methods = list(
    initialize = function(boundary) {
      boundary <<- boundary
    },
    
    randomize = function(n) {
      require(sp)
      spcoord <- spsample(boundary, n, type="random", iter=10)
      if (is.null(spcoord)) return(NULL)
      return(spcoord)
    },
    
    plotLocations = function(locations) {
      require(sp)
      plot(boundary)
      plot(locations, add=T, col="blue")
    }
  )
)

ClusteredInitialPopulation <- setRefClass(
  Class = "ClusteredInitialPopulation",
  contains = "InitialPopulation",
  fields = list(
    studyArea = "StudyArea",
    habitatWeights = "HabitatWeights",
    weights = "numeric",
    mesh = "inla.mesh",
    spde = "inla.spde2",
    Q = "Matrix",
    samples = "numeric"
  ),
  
  methods = list(
    initialize = function(studyArea, habitatWeights=HabitatWeights(), max.edge=5000, mu=-22, sigma=1, range=100*1e3, fun=exp, seed=1L, ...) {
      callSuper(...)
      studyArea <<- studyArea
      habitatWeights <<- habitatWeights
      
      sampleMaternRandomField(range=range, mu=mu, sigma=sigma, seed=seed, max.edge=max.edge, fun=fun)
      if (class(habitatWeights) == "HabitatWeights") weights <<- 1
      else setHabitatWeightsForField()      
    },
    
    sampleMaternRandomField = function(range=100e3, mu=-22, sigma=1, seed=1L, max.edge=5000, fun=exp) {
      message("Sampling from MRF with range = ", range, ", sigma = ", sigma, ", mu = ", mu, ", seed = ", seed, "...")
      fun <- match.fun(fun)
      mesh <<- inla.mesh.create(boundary=inla.sp2segment(studyArea$boundary), refine=list(max.edge=max.edge))
      message("Nodes for randomizing individual locations = ", mesh$n)
      kappa <- sqrt(8) / range
      spde <<- inla.spde2.matern(mesh, alpha=2)
      theta <- c(-0.5 * log(4 * pi * sigma^2 * kappa^2), log(kappa))
      #message("theta1 = ", theta[1], ", theta2 = ", theta[2])
      Q <<- inla.spde2.precision(spde, theta)
      samples <<- fun(as.vector(inla.qsample(mu=rep(mu, times=nrow(Q)), Q=Q, seed=seed)))
      #message("sample mean = ", mean(samples), " sd = ", sd(samples))
    },
    
    setHabitatWeightsForField = function() {
      message("Extracting habitat types...")
      locations <- mesh$loc[,1:2,drop=F]
      habitatTypes <- extract(studyArea$habitat, locations)
      weights <<- habitatWeights$getWeights(habitatTypes)
    },
    
    randomize = function(n) {
      require(sp)
      index <- sample(1:mesh$n, size=n, replace=T, prob=samples * weights)
      locations <- mesh$loc[index,]
      return(SpatialPoints(locations[,1:2,drop=F], proj4string=studyArea$proj4string))      
    },
    
    plotMesh = function() {
      plot(mesh)
    },
    
    plotLocations = function(locations) {
      require(lattice)
      require(grid)
      
      proj <- inla.mesh.projector(mesh, dims=dim(studyArea$habitat)[1:2] / studyArea$plotScale)
      z <- inla.mesh.project(proj, field=samples * weights)
      
      panel <- if (missing(locations)) panel.levelplot
      else {
        xy <- coordinates(locations)
        function(...) {
          panel.levelplot(...)
          grid.points(xy[,1], xy[,2], pch='.')
        }
      }
      
      p <- levelplot(row.values=proj$x, column.values=proj$y, x=z, panel=panel)
      print(p)
    }
  )
)
