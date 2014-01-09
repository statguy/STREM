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
    habitatWeigts = "HabitatWeights",
    weights = "numeric",
    mesh = "inla.mesh",
    spde = "inla.spde2",
    Q = "Matrix",
    sample = "numeric",
    probabilities = "numeric"
  ),
  
  methods = list(
    initialize = function(studyArea, habitatWeights, max.edge=5000, range=100*1e3, sigma=1, seed=0L, ...) {
      callSuper(...)
      studyArea <<- studyArea
      habitatWeights <<- habitatWeights
      
      randomizeMaternField(max.edge=max.edge, range=range, sigma=sigma, seed=seed)
      setHabitatWeightsForField()      
    },
    
    randomizeMaternField = function(max.edge=5000, range=100*1e3, sigma=1, seed=0L) {
      mesh <<- inla.mesh.create(boundary=inla.sp2segment(studyArea$boundary), refine=list(max.edge=max.edge))
      message("Nodes for randomizing individual locations = ", mesh$n)
      message("Sampling...")
      sigma0 <- sigma
      range0 <- range
      kappa0 <- sqrt(8)/range0
      tau0 <- 1/(sqrt(4*pi)*kappa0*sigma0)
      spde <<- inla.spde2.matern(mesh,
                                 B.tau=cbind(log(tau0),1,0),
                                 B.kappa=cbind(log(kappa0),0,1),
                                 theta.prior.mean=c(0,0),
                                 theta.prior.prec=1)
      Q <<- inla.spde2.precision(spde, theta=c(0,0))
      sample <<- exp(as.vector(inla.qsample(n=1, Q=Q, seed=seed)))
      probabilities <<- (sample-min(sample))/sum(sample-min(sample))
    },
    
    setHabitatWeightsForField = function() {
      message("Extracting habitat types...")
      locations <- mesh$loc[,1:2,drop=F]
      habitatTypes <- extract(studyArea$habitat, locations)
      weights <<- habitatWeights$getWeights(habitatTypes)
    },
    
    randomize = function(n) {
      require(sp)
      index <- sample(1:mesh$n, size=n, replace=T, prob=probabilities * weights)
      locations <- mesh$loc[index,]
      return(SpatialPoints(locations[,1:2,drop=F], proj4string=studyArea$proj4string))      
    },
    
    plotMesh = function() {
      plot(mesh)
    },
    
    plotLocations = function(locations) {
      require(lattice)

      proj <- inla.mesh.projector(mesh, dims=dim(studyArea$habitat)[1:2]/100)
      z <- inla.mesh.project(proj, field=probabilities * weights)
      
      p <- if (missing(locations)) {
        levelplot(row.values=proj$x, column.values=proj$y, x=z)
      }
      else {
        require(sp)
        require(grid)
        xy <- coordinates(locations)
        levelplot(row.values=proj$x, column.values=proj$y, x=z)
        trellis.focus("panel", 1, 1, highlight=FALSE)
        lpoints(xy[,1], xy[,2], pch='.')
        trellis.unfocus()
      }
      print(p)
    }
  )
)
