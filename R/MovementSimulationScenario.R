MovementSimulationScenario <- setRefClass(
  Class = "MovementSimulationScenario",
  #contains = "SimulationStudy",
  fields = list(
    study = "SimulationStudy",
    
    years = "integer",
    days = "integer",
    stepIntervalHours = "numeric",
    stepSpeedScale = "numeric",
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
    initialize = function(...) {
      distanceScale <<- 1e3
      loadHabitatRasterInMemory <<- FALSE
      callSuper(...)
      return(.self)
    },
    
    newInstance = function() {
      nSteps.tmp <- 24 * days / stepIntervalHours
      message("Number of steps = ", nSteps.tmp, ", steps per day = ", 24 / stepIntervalHours)
      if (nSteps.tmp %% 1 != 0) stop("Number of steps must be integer.")
      nSteps <<- as.integer(nSteps.tmp)
      
      #studyArea <<- if (isTest) TestStudyArea$new(context=context)$newInstance()
      #else FinlandStudyArea$new(context=context)$newInstance()
      return(.self)
    },
    
    finalize = function() {
    },
        
    randomizeVector = function(locations) {
      library(sp)
      library(raster)
      
      point <- SpatialPoints(locations[,,drop=F], proj4string=study$studyArea$proj4string)
      outsideBoundary <- is.na(over(point, study$studyArea$boundary))
      if (class(outsideBoundary) == "matrix") outsideBoundary <- outsideBoundary[,1]

      if (inherits(habitatWeights, "uninitializedField")) {
        if (all(outsideBoundary == TRUE)) return(NULL)
        return(list(index=1, coords=locations[1,,drop=F]))
      }
      else {
        habitatTypes <- extract(study$studyArea$habitat, locations)
        w <- habitatWeights$getWeights(habitatTypes)
        w[outsideBoundary] <- 0
        
        if (all(w==0)) return(NULL)
        k <- sample(1:nrow(locations), 1, prob=w)
        return(list(index=k, coords=locations[k,,drop=F]))
      }
    },
        
    randomizeBCRWTrack = function(initialLocation, initialAngle, isFirst, nProposal) {
      library(CircStats)
      library(sp)
      
      maxTry <- 100
            
      coords <- matrix(NA, nrow=nSteps + 1, ncol=2)
      coords[1,] <- initialLocation
      angles <- numeric(nSteps + 1)
      angles[1] <- initialAngle

      step <- 2
      while (TRUE) {
        for (j in 1:maxTry) {
          # Correlated random walk
          newAngles <- rwrpnorm(nProposal, angles[step-1], CRWCorrelation)
          newDistances <- rweibull(nProposal, shape=2, scale=stepSpeedScale * stepIntervalHours * distanceScale)
          
          if (length(homeRangeRadius) != 0) {
            if (step > 2 & euclidean(coords[1,,drop=F], coords[step-1,,drop=F]) > homeRangeRadius) {
              # Biased correlated random walk
              xy <- matrix(coords[1,,drop=F], ncol=2, nrow=nProposal, byrow=T) - matrix(coords[step-1,,drop=F], ncol=2, nrow=nProposal, byrow=T)
              angleBias <- atan2(xy[,2], xy[,1])
              newAngles <- Arg(BCRWCorrelationBiasTradeoff * exp(complex(imaginary=newAngles)) + (1-BCRWCorrelationBiasTradeoff) * exp(complex(imaginary=angleBias)))
            }
          }
          
          proposedVectors <- getVector(coords[step-1,,drop=F], newDistances, newAngles)
          acceptedVectors <- randomizeVector(proposedVectors)
          
          if (!is.null(acceptedVectors)) {
            coords[step,] <- acceptedVectors$coords
            angles[step] <- newAngles[acceptedVectors$index]
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
          track <- SpatialLines(list(Lines(list(Line(coords[1:step-1,])), ID=1)), proj4string=study$studyArea$proj4string)
          plot(track, col="blue")
          plot(study$studyArea$boundary, add=T)
          points(proposedVectors, col="red", pch="+")
          points(acceptedVectors$coords, col="green", pch="+")
          stop("Boundary reflection failed.")
        }
        
        if (step == nSteps + 1) break
        step <- step + 1
      }
      
      index <- if (isFirst) 1:nSteps else 1:nSteps+1
      stepDays <- rep(1:days, each=24 / stepIntervalHours)
      stepHours <- rep(seq(0, 24-stepIntervalHours, by=stepIntervalHours), days)
      stepMinutes <- (stepHours - trunc(stepHours)) * 60
      stepSeconds <- (stepMinutes - trunc(stepMinutes)) * 60
      return(data.frame(x=coords[index,1], y=coords[index,2], angle=angles[index], day=stepDays, hour=trunc(stepHours), minute=trunc(stepMinutes), second=round(stepSeconds)))
    },
    
    # Combined birth-death process:
    # 0 born = individual dies/emigrates
    # 1 born = individual survives to the next year
    # 2 born = 1 parent survives + 1 offspring/immigration
    # etc.
    randomizeBirthDeath = function(param=list(mean=1, sd=1.1)) {
      bdRate <- rlnorm(n=1, meanlog=log(param$mean), sdlog=log(param$sd))
      nTransform <- rpois(length(agents), bdRate)
      
      nBorn <- sum(nTransform[nTransform > 1] - 1)
      nSurvive <- sum(nTransform[nTransform == 1]) + length(nTransform[nTransform > 1])
      nDie <- length(nTransform[nTransform == 0])
      
      message("agents before = ", length(agents), " -> born = ", nBorn, ", survive = ", nSurvive, ", die = ", nDie, " -> agents after = ", nSurvive + nBorn)
      
      survivedIndex <- nTransform > 0      
      bornIndex <- nTransform > 1
      x <- rep(which(bornIndex), nTransform[bornIndex]-1)
      survivedBornIndex <- c(which(survivedIndex), x)
      
      agents <<- agents[survivedIndex]
      if (nBorn > 0) agents <<- c(agents, newAgentId:(newAgentId + nBorn - 1))
      newAgentId <<- as.integer(newAgentId + nBorn)
      
      return(survivedBornIndex)
    },
    
    randomizeBCRWTracks = function() {
      library(plyr)
      library(maptools)
      
      initialLocations <- initialPopulation$randomize(nAgents)
      habitatTypes <- extract(study$studyArea$habitat, initialLocations)
      if (any(is.na(habitatTypes))) stop("Invalid initial coordinates.")
      
      nProposal <- if (inherits(habitatWeights, "uninitializedField")) 1 else 10
      message("Number of proposals = ", nProposal)
      
      agents <<- 1:nAgents
      newAgentId <<- as.integer(nAgents + 1)
      
      tracks <- data.frame()
      initialLocations <- coordinates(initialLocations)
      initialAngles <- runif(nAgents, 0, 2*pi)
      isFirst <- TRUE
      
      nAgentsCurrent <- nAgents
      
      for (year in 1:years) {
        if (nAgentsCurrent > 0) {
          track <- ldply(1:nAgentsCurrent,
            function(agentIndex, initialLocations, initialAngles, isFirst, nProposal) {
              library(sp)
              
              ## TODO: UNTESTED CODE
              if (length(habitatWeights) != 0 & loadHabitatRasterInMemory) {
                study$studyArea$readRasterIntoMemory()                
              }
              
              message("Agent (", agents[agentIndex], ") = ", agentIndex, " / ", nAgentsCurrent, ", Year = ", year,  " / ", years, "...")
              track <- randomizeBCRWTrack(initialLocation=initialLocations[agentIndex,,drop=F], initialAngle=initialAngles[agentIndex], isFirst=isFirst, nProposal=nProposal)
              track$agent <- agents[agentIndex]
              return(track)
            },
            initialLocations=initialLocations, initialAngles=initialAngles, isFirst=isFirst, nProposal=nProposal, .parallel=runParallel)
          
          track$year <- year
          tracks <- rbind(tracks, track)
        }
        
        if (year < years) {
          survivedBornLastStepIndex <- randomizeBirthDeath() * nSteps
          initialLocations <- as.matrix(track[survivedBornLastStepIndex, c("x","y"), drop=F])
          initialAngles <- track[survivedBornLastStepIndex, c("angle")]
          isFirst <- FALSE
          nAgentsCurrent <- length(agents)
        }
      }

      return(tracks)
    },
    
    simulate = function(restartIteration=1, save=FALSE) {
      stopifnot(restartIteration <= nIterations)
      
      simulatedTracks <- SimulatedTracksCollection$new(study=study)
      
      for (i in restartIteration:nIterations) {
        message("Iteration ", i, " of ", nIterations, "...")
        tracksDF <- randomizeBCRWTracks()
        date <- as.POSIXct(strptime(paste(2000+tracksDF$year, tracksDF$day, tracksDF$hour, tracksDF$minute, tracksDF$second), format="%Y %j %H %M %S"))
        tracks <- SimulatedTracks$new(study=study, preprocessData=save, xy=tracksDF[,c("x","y")], id=tracksDF$agent, date=date, iteration=i)
        simulatedTracks$add(tracks)
      }
      
      return(invisible(simulatedTracks))
    }        
  )
)

# 6 minutes step length for distance fix
MovementSimulationScenarioIntensive <- setRefClass(
  Class = "MovementSimulationScenarioIntensive",
  contains = "MovementSimulationScenario",
  methods = list(
    initialize = function(nAgents=as.integer(2), nIterations=as.integer(1), years=as.integer(2), days=as.integer(10), stepIntervalHours=0.1, runParallel=T, isTest=F, ...) {
      callSuper(years=years, nAgents=nAgents, nIterations=nIterations, days=days, stepIntervalHours=stepIntervalHours, stepSpeedScale=0.5, CRWCorrelation=0.8, runParallel=runParallel, ...)
      return(.self)
    },
    
    newInstance = function(context, response="Intensive", isTest=F) {
      callSuper()
      study <<- SimulationStudy$new(response=response)$newInstance(context=context, isTest=isTest)
      initialPopulation <<- RandomInitialPopulation$new(studyArea=study$studyArea)
      return(.self)
    }
  )
)

if (F) {

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
      initialPopulation <<- RandomInitialPopulation$new(studyArea=study$studyArea)
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
    #  library(plyr)
    #  x <- ddply(trackReplicates, .(year, individual, iteration), function(x, i) x[i,], i=retainDaysIndex)
    #  return(x)
    #}
  )
)

# TODO: scenario C

# Biased correlated random walk in a homogenous landscape, spatial variation in initial locations
MovementSimulationScenarioD <- setRefClass(
  Class = "MovementSimulationScenarioD",
  contains = "MovementSimulationScenarioA",
  fields = list(
    
  ),
  methods = list(
    initialize = function(isTest=FALSE, response="D", ...) {
      callSuper(response=response, isTest=isTest, ...)
      initialPopulation <<- ClusteredInitialPopulation$new(studyArea=study$studyArea, habitatWeights=CORINEHabitatWeights())
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
      habitatWeights <<- CORINEHabitatWeights$new(list(urban=0.1, agriculture=0.1, forestland=1, peatland=0.5, water=0.05))
      return(.self)
    }
  )
)

}