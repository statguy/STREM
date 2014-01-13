MovementSimulationScenario <- setRefClass(
  Class = "MovementSimulationScenario",
  fields = list(
    context = "ANY",
    response = "character",
    years = "integer",
    days = "integer",
    stepIntervalHours = "numeric",
    stepSpeedScale = "numeric",
    isTest = "logical",
    studyArea = "StudyArea",
    initialPopulation = "InitialPopulation",
    habitatWeights="ANY",
    nAgents = "integer",
    nIterations = "integer",
    CRWCorrelation = "numeric",
    BCRWCorrelationBiasTradeoff = "numeric",
    homeRangeRadius = "numeric",
    distanceScale = "numeric",
    loadHabitatRasterInMemory = "logical",
    runParallel = "logical"
  ),
  methods = list(
    initialize = function(isTest=FALSE, ...) {
      distanceScale <<- 1e3
      loadHabitatRasterInMemory <<- FALSE
      callSuper(...)
      isTest <<- isTest

      # TODO: start cluster
      
      return(.self)
    },
    
    newInstance = function() {
      studyArea <<- if (isTest) TestStudyArea(context=context)$newInstance()
      else FinlandStudyArea(context=context)$newInstance()
      return(.self)
    },
    
    finalize = function() {
      # TODO: stop cluster 
    },
        
    euclidean = function(c1, c2) sqrt( (c1[,1] - c2[,1])^2 + (c1[,2] - c2[,2])^2 ),
    
    newVector = function(coords, distance, angle) {
      return(cbind(coords[,1] + distance * cos(angle), coords[,2] + distance * sin(angle)))
    },
    
    randomizeVector = function(locations)
    {
      require(sp)
      require(raster)
      
      point <- SpatialPoints(locations[,,drop=F], proj4string=studyArea$proj4string)
      outsideBoundary <- is.na(over(point, studyArea$boundary))
      if (class(outsideBoundary) == "matrix") outsideBoundary <- outsideBoundary[,1]

      if (class(habitatWeights) == "uninitializedField") {
        if (all(outsideBoundary == TRUE)) return(NULL)
        return(list(index=1, coords=locations[1,,drop=F]))
      }
      else {
        habitatTypes <- extract(studyArea$habitat, locations)
        w <- habitatWeights$getWeights(habitatTypes)
        w[outsideBoundary] <- 0
        
        if (all(w==0)) return(NULL)
        k <- sample(1:nrow(locations), 1, prob=w)
        return(list(index=k, coords=locations[k,,drop=F]))
      }
    },

    randomizeBCRWTrack = function(initialLocation, nProposal) {
      require(CircStats)
      require(sp)
      
      maxTry <- 100
      nSteps <- 24 * days / stepIntervalHours
      
      coords <- matrix(NA, nrow=nSteps, ncol=2)
      coords[1,] <- initialLocation
      prevAngle <- runif(1, 0, 2*pi)    
                
      for (step in 2:nSteps) {
        for (j in 1:maxTry) {
          # Correlated random walk
          newAngles <- rwrpnorm(nProposal, prevAngle, CRWCorrelation)
          newDistances <- rweibull(nProposal, shape=2, scale=stepSpeedScale * stepIntervalHours * distanceScale)
          
          if (length(homeRangeRadius) != 0) {
            if (step > 2 & euclidean(coords[1,,drop=F], coords[step-1,,drop=F]) > homeRangeRadius) {
              # Biased correlated random walk
              xy <- matrix(coords[1,,drop=F], ncol=2, nrow=nProposal, byrow=T) - matrix(coords[step-1,,drop=F], ncol=2, nrow=nProposal, byrow=T)
              angleBias <- atan2(xy[,2], xy[,1])
              newAngles <- Arg(BCRWCorrelationBiasTradeoff * exp(complex(imaginary=newAngles)) + (1-BCRWCorrelationBiasTradeoff) * exp(complex(imaginary=angleBias)))
            }
          }
          
          v <- newVector(coords[step-1,,drop=F], newDistances, newAngles)
          indexVector <- randomizeVector(v)
          
          if (!is.null(indexVector)) {
            prevAngle <- newAngles[indexVector$index]
            coords[step,] <- indexVector$coords
            break
          }
          else {
            # Randomize angle at the boundary, so we don't get stuck
            prevAngle <- runif(1, 0, 2*pi)
            next
          }
        }
        
        if (j == maxTry) {
          track <- SpatialLines(list(Lines(list(Line(coords[1:step-1,])))), proj4string=studyArea$proj4string)
          plot(track, col="blue")
          plot(boundary, add=T)
          points(v, col="red", pch="+")
          points(k, col="green", pch="+")
          stop("Boundary reflection failed. This should not happen.")
        }
      }
      
      return(data.frame(x=coords[,1], y=coords[,2], year = rep(1:years, each=nSteps)))
    },
    
    randomizeBCRWTracks = function() {
      require(plyr)
      
      initialLocations <- initialPopulation$randomize(nAgents)
      habitatTypes <- extract(studyArea$habitat, initialLocations)
      if (any(is.na(habitatTypes))) stop("Invalid initial coordinates.")
      
      nProposal <- if (length(habitatWeights) == 0) 1 else 10
      message("Number of proposals = ", nProposal)
      
      nSteps <- 24 * days / stepIntervalHours
      message("Number of steps = ", nSteps, ", steps per day = ", 24 / stepIntervalHours)
      if (nSteps %% 1 != 0) stop("Number of steps must be integer.")

      initialLocations <- coordinates(initialLocations)
      tracks <- ldply(1:nAgents,
        function(i, initialLocations, nProposal) {
          require(sp)
          message("Agent = ", i, " / ", nrow(initialLocations), "...")
          
          ## TODO: UNTESTED CODE
          if (length(habitatWeights) != 0 & loadHabitatRasterInMemory) {
            require(raster)
            require(rgdal)
            if (!inMemory(studyArea$habitat)) {
              message("Reading habitat raster to memory...")
              studyArea$habitat <<- readAll(studyArea$habitat)
            }
          }
          
          track <- randomizeBCRWTrack(initialLocation=initialLocations[i,,drop=F], nProposal=nProposal)
          track$individual <- i
                    
          return(track)
        },
        initialLocations=initialLocations, nProposal=nProposal, .parallel=runParallel)
      
      return(tracks)
    },

    simulate = function() {
      trackReplicates <- data.frame()
      for (i in 1:nIterations) {
        message("Iteration ", i, " of ", nIterations, "...")
        tracks <- randomizeBCRWTracks()
        tracks$iteration <- i
        trackReplicates <- rbind(trackReplicates, tracks)
      }
      return(trackReplicates)
    },
    
    plotTracks = function(tracks) {
      require(sp)
      require(plyr)
      plot(studyArea$boundary)
      ddply(tracks, .(individual, year, iteration), function(tracks) {
        lines(tracks$x, tracks$y, col=tracks$individual)        
      })
      invisible()
    },
    
    getTracksFile = function() {
      return(context$getFileName(context$resultDataDirectory, name="SimulatedTracks", response=response, region=studyArea$region))
    },
    
    saveTracks = function(tracks) {
      save(tracks, file=getTracksFile())
    },
    
    loadTracks = function() {
      load(getTracksFile())
      return(tracks)
    }    
  )
)

# Correlated random walk in a homogeneous landscape, random initial locations
MovementSimulationScenarioA <- setRefClass(
  Class = "MovementSimulationScenarioA",
  contains = "MovementSimulationScenario",
  methods = list(
    initialize = function(isTest=FALSE, response="A", ...) {
      callSuper(response=response, isTest=isTest, ...)
      return(.self)
    },

    newInstance = function() {
      callSuper()
      initialPopulation <<- RandomInitialPopulation(boundary=studyArea$boundary)
      return(.self)
    }
  )
)

# Biased correlated random walk in a homogenous landscape, random initial locations
MovementSimulationScenarioB <- setRefClass(
  Class = "MovementSimulationScenarioB",
  contains = "MovementSimulationScenarioA",
  fields = list(

  ),
  methods = list(
    initialize = function(isTest=FALSE, response="B", ...) {
      callSuper(response=response, isTest=isTest, ...)
    }
    
    #discardMovements = function(tracks, retainDaysIndex) {
    #  require(plyr)
    #  x <- ddply(trackReplicates, .(year, individual, iteration), function(x, i) x[i,], i=retainDaysIndex)
    #  return(x)
    #}
  )
)

# Biased correlated random walk in a homogenous landscape, spatial variation in initial locations
MovementSimulationScenarioD <- setRefClass(
  Class = "MovementSimulationScenarioD",
  contains = "MovementSimulationScenarioA",
  fields = list(
    
  ),
  methods = list(
    initialize = function(isTest=FALSE, response="D", ...) {
      callSuper(response=response, isTest=isTest, ...)
      initialPopulation <<- ClusteredInitialPopulation(boundary=studyArea$boundary, habitatWeights=CORINEHabitatWeights())
    }
  )
)

# Correlated random walk in a heterogenous landscape, random initial locations
MovementSimulationScenarioE <- setRefClass(
  Class = "MovementSimulationScenarioE",
  contains = "MovementSimulationScenarioA",
  methods = list(
    initialize = function(isTest=FALSE, response="E", ...) {
      callSuper(response=response, isTest=isTest, ...)
      return(.self)
    },
    
    newInstance = function() {
      callSuper()
      habitatWeights <<- CORINEHabitatWeights(list(urban=0.1, agriculture=0.1, forestland=1, peatland=0.5, water=0.05))
      return(.self)
    }
  )
)

