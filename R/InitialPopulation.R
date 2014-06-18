library(INLA)
setOldClass("inla.mesh")
setOldClass("inla.spde2")

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
      return(invisible(.self))
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
      library(sp)
      locations <<- spsample(studyArea$boundary, n, type="random", iter=10)
      return(locations)
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
    initialize = function(studyArea, habitatWeights=HabitatWeights(), max.edge=5000, sigma=1, range=100*1e3, fun=exp, seed=1L, ...) {
      callSuper(...)
      studyArea <<- studyArea
      habitatWeights <<- habitatWeights
      
      sampleMaternRandomField(range=range, sigma=sigma, seed=seed, max.edge=max.edge, fun=fun)
      if (inherits(habitatWeights, "HabitatWeights")) weights <<- 1
      else setHabitatWeightsForField()
    },
    
    sampleMaternRandomField = function(range=100e3, sigma=1, seed=1L, max.edge=5000, fun=exp) {
      library(INLA)
      message("Sampling from GMRF with range = ", range, ", sigma = ", sigma, ", seed = ", seed, "...")
      fun <- match.fun(fun)
      mesh <<- inla.mesh.create(boundary=inla.sp2segment(studyArea$boundary), refine=list(max.edge=max.edge))
      message("Nodes for randomizing individual locations = ", mesh$n)
      kappa <- sqrt(8) / range
      spde <<- inla.spde2.matern(mesh, alpha=2)
      theta <- c(-0.5 * log(4 * pi * sigma^2 * kappa^2), log(kappa))
      Q <<- inla.spde2.precision(spde, theta)
      samples <<- fun(as.vector(inla.qsample(Q=Q, seed=seed)))
      #message("sample mean = ", mean(samples), " sd = ", sd(samples))
      return(invisible(.self))
    },
    
    setHabitatWeightsForField = function() {
      library(raster)
      message("Extracting habitat types...")
      meshLocations <- mesh$loc[,1:2,drop=F]
      habitatTypes <- extract(studyArea$habitat, meshLocations)
      weights <<- habitatWeights$getWeights(habitatTypes)
      return(invisible(.self))
    },
    
    randomize = function(n) {
      library(sp)
      index <- sample(1:mesh$n, size=n, replace=T, prob=samples * weights)
      xy <- mesh$loc[index,]
      locations <<- SpatialPoints(xy[,1:2,drop=F], proj4string=studyArea$proj4string)
      return(locations)
    },
    
    plotMesh = function() {
      library(INLA)
      plot(mesh)
      return(invisible(.self))
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
      
      return(invisible(.self))
    }
  )
)
