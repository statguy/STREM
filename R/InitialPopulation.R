library(INLA)
setOldClass("inla.mesh")
setOldClass("inla.spde2")
library(sp)

InitialPopulation <- setRefClass(
  Class = "InitialPopulation",
  fields = list(
    studyArea = "StudyArea",
    locations = "SpatialPoints"
  ),
  methods = list(
    randomize = function(n) {
      stop("randomize() undefined.")
    },
    
    plotLocations = function() {
      library(sp)
      plot(locations, col="blue")
      plot(studyArea$boundary, add=T)
    }
  )
)

RandomInitialPopulation <- setRefClass(
  Class = "RandomInitialPopulation",
  contains = "InitialPopulation",
  fields = list(
  ),
  methods = list(
    randomize = function(n) {
      locations <<- spsample(studyArea$boundary, n, type="random", iter=10)
      invisible(locations)
    }
  )
)

ClusteredInitialPopulation <- setRefClass(
  Class = "ClusteredInitialPopulation",
  contains = "InitialPopulation",
  fields = list(
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
      if (inherits(habitatWeights, "HabitatWeights")) weights <<- 1
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
      meshLocations <- mesh$loc[,1:2,drop=F]
      habitatTypes <- extract(studyArea$habitat, mehsLocations)
      weights <<- habitatWeights$getWeights(habitatTypes)
    },
    
    randomize = function(n) {
      index <- sample(1:mesh$n, size=n, replace=T, prob=samples * weights)
      xy <- mesh$loc[index,]
      locations <<- SpatialPoints(xy[,1:2,drop=F], proj4string=studyArea$proj4string)
      invisible(locations)
    },
    
    plotMesh = function() {
      plot(mesh)
    },
    
    plotLocations = function() {
      library(lattice)
      library(grid)
      
      proj <- inla.mesh.projector(mesh, dims=dim(studyArea$habitat)[1:2] / studyArea$plotScale)
      z <- inla.mesh.project(proj, field=samples * weights)
      
      panel <- if (inherits(locations, "uninitializedField")) panel.levelplot
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
