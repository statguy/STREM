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
    runParallel = "logical",
    
    nSteps = "integer",
    birthRate = "numeric",
    deathRate = "numeric",
    agents = "integer",
    newAgentId = "integer"
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
      nSteps.tmp <- 24 * days / stepIntervalHours
      message("Number of steps = ", nSteps.tmp, ", steps per day = ", 24 / stepIntervalHours)
      if (nSteps.tmp %% 1 != 0) stop("Number of steps must be integer.")
      nSteps <<- as.integer(nSteps.tmp)
      
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
    
    randomizeVector = function(locations) {
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
      
      coords <- matrix(NA, nrow=nSteps, ncol=2)
      coords[1,] <- initialLocation
      prevAngle <- runif(1, 0, 2*pi)    
      
      step <- 2
      while (TRUE) {
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
          
          proposedVectors <- newVector(coords[step-1,,drop=F], newDistances, newAngles)
          acceptedVectors <- randomizeVector(proposedVectors)
          
          if (!is.null(acceptedVectors)) {
            prevAngle <- newAngles[acceptedVectors$index]
            coords[step,] <- acceptedVectors$coords
            break
          }
          else {
            # Vector is pointing outside the boundary or there is only unknown habitat.
            # Go back one step and try again, so we don't get stuck.
            step <- step - 1
            if (step < 2) step <- 2
            next
          }
        }
        
        if (j == maxTry) {
          track <- SpatialLines(list(Lines(list(Line(coords[1:step-1,])), ID=1)), proj4string=studyArea$proj4string)
          plot(track, col="blue")
          plot(studyArea$boundary, add=T)
          points(proposedVectors, col="red", pch="+")
          points(acceptedVectors$coords, col="green", pch="+")
          stop("Boundary reflection failed.")
        }

        if (step == nSteps) break
        step <- step + 1
      }
      
      return(data.frame(x=coords[,1], y=coords[,2]))
    },
    
    randomizeBirths = function() {
      nParents <- 0
      for (agent in agents) {  
        nBorn <- rpois(1, birthRate)
        if (nBorn>0) {
          nParents <- nParents + 1
          agents <<- c(agents, newAgentId:(newAgentId + nBorn - 1))
          newAgentId <<- as.integer(newAgentId + nBorn)
        }
      }

      message("Born = ", length(agents) - nAgents, " for = ", nParents, " out of ", nAgents)
      
      nAgents <<- length(agents)
      if (nAgents > 2000) stop("Too many agents to simulate.")
    },
    
    randomizeDeaths = function() {
      deathIndex <- runif(nAgents) <= 1 - exp(-deathRate * nAgents)
      agents <<- agents[!deathIndex]

      message("Died = ", nAgents - length(agents), " out of ", nAgents)
      
      nAgents <<- length(agents)
      if (nAgents == 0) stop("Agents extint.")
    },
    
    randomizeBCRWTracks = function() {
      require(plyr)
      require(maptools)
      
      initialLocations <- initialPopulation$randomize(nAgents)
      habitatTypes <- extract(studyArea$habitat, initialLocations)
      if (any(is.na(habitatTypes))) stop("Invalid initial coordinates.")
      
      nProposal <- if (class(habitatWeights) == "uninitializedField") 1 else 10
      message("Number of proposals = ", nProposal)
      
      agents <<- 1:nAgents
      newAgentId <<- as.integer(nAgents + 1)
      
      tracks <- data.frame()
      initialLocations <- coordinates(initialLocations)
      
      for (year in 1:years) {
        track <- ldply(agents,
          function(agentId, initialLocations, nProposal) {
            require(sp)

            ## TODO: UNTESTED CODE
            if (length(habitatWeights) != 0 & loadHabitatRasterInMemory) {
              require(raster)
              require(rgdal)
              if (!inMemory(studyArea$habitat)) {
                message("Reading habitat raster to memory...")
                studyArea$habitat <<- readAll(studyArea$habitat)
              }
            }
            
            message("Agent = ", agentId, " / ", nAgents, ", Year = ", year,  " / ", years, "...")
            track <- randomizeBCRWTrack(initialLocation=initialLocations[agentId,,drop=F], nProposal=nProposal)
            track$agent <- agentId
            return(track)
          },
          initialLocations=initialLocations, nProposal=nProposal, .parallel=runParallel)

        track$year <- year
        tracks <- rbind(tracks, track)
        initialLocations <- as.matrix(track[seq(nSteps, nAgents*nSteps, by=nSteps), c("x","y"), drop=F])
        
        if (year < years) {
          #randomizeBirths()
          #randomizeDeaths()
        }
      }

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
      
      spTracks <- dlply(trackReplicates, .(agent, year, iteration), function(x) {
        id <- paste(x[1, c("agent","year","iteration")], collapse=".")
        return(Lines(list(Line(x[,c("x","y")])), ID=id))
      }, .parallel=runParallel)
      spData <- ddply(trackReplicates, .(agent, year, iteration), function(x) {
        return(x[1, c("agent","year","iteration"), drop=F])
      }, .parallel=runParallel)
      spdfTracks <- SpatialLinesDataFrame(SpatialLines(spTracks, proj4string=studyArea$proj4string),
                                          spData,
                                          match.ID=FALSE)      
      return(spdfTracks)
    },
    
    plotTracks = function(tracks) {
      require(sp)
      plot(tracks, col=tracks@data$agent)
      plot(studyArea$boundary, add=T, border="lightgray")
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

# 6 minutes step length for distance fix
MovementSimulationScenarioIntensive <- setRefClass(
  Class = "MovementSimulationScenarioIntensive",
  contains = "MovementSimulationScenario",
  methods = list(
    initialize = function(response="Intensive", nAgents=as.integer(100), nIterations=as.integer(1), runParallel=T, ...) {
      callSuper(response=response, isTest=F, years=as.integer(1), nAgents=nAgents, nIterations=nIterations, days=as.integer(60), stepIntervalHours=0.1, stepSpeedScale=0.5, CRWCorrelation=0.8, runParallel=runParallel, ...)
      return(.self)
    },
    
    newInstance = function() {
      callSuper()
      initialPopulation <<- RandomInitialPopulation(boundary=studyArea$boundary)
      return(.self)
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
