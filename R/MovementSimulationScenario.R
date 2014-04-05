MovementSimulationScenario <- setRefClass(
  Class = "MovementSimulationScenario",
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
    BCRWCorrelationBiasTradeoff = "ANY",
    homeRangeRadius = "ANY",
    distanceScale = "numeric",
    loadHabitatRasterInMemory = "logical",
    runParallel = "logical",
    
    nSteps = "integer",
    agents = "integer",
    newAgentId = "integer",
    nProposal = "integer"
  ),
  methods = list(
    initialize = function(...) {
      distanceScale <<- 1e3
      loadHabitatRasterInMemory <<- FALSE
      callSuper(...)
      return(invisible(.self))
    },
    
    newInstance = function() {
      BCRWCorrelationBiasTradeoff <<- rep(NA, nAgents)
      homeRangeRadius <<- rep(NA, nAgents)
      nSteps.tmp <- 24 * days / stepIntervalHours
      message("Number of steps = ", nSteps.tmp, ", steps per day = ", 24 / stepIntervalHours)
      if (nSteps.tmp %% 1 != 0) stop("Number of steps must be integer.")
      nSteps <<- as.integer(nSteps.tmp)
      return(invisible(.self))
    },
    
    randomizeDistance = function(n) {
      rweibull(n, shape=2, scale=stepSpeedScale * stepIntervalHours * distanceScale)
    },
    
    randomizeVector = function(locations) {
      library(sp)
      library(raster)
      
      if (nrow(locations) < 1 | ncol(locations) != 2 | any(is.nan(locations))) {
        print(locations)
        stop("Invalid locations argument.")
      }  
      
      point <- SpatialPoints(locations[,,drop=F], proj4string=study$studyArea$boundary@proj4string)
      
      if (inherits(point, "numeric")) {
        print(locations)
        print(point)
        stop("Invalid locations argument.")
      }
      
      outsideBoundary <- is.na(over(point, study$studyArea$boundary))
      if (class(outsideBoundary) == "matrix") outsideBoundary <- outsideBoundary[,1]
      
      if (inherits(habitatWeights, "uninitializedField") | is.null(habitatWeights)) {
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
    
    randomizeBCRWTrack = function(initialLocation, initialAngle, isFirst, CRWCorrelation, BCRWCorrelationBiasTradeoff, homeRangeRadius) {
      library(CircStats)
      library(sp)
      
      maxTry <- 10000
      
      coords <- matrix(NA, nrow=nSteps + 1, ncol=2)
      coords[1,] <- initialLocation
      angles <- numeric(nSteps + 1)
      angles[1] <- initialAngle
      
      step <- 2
      
      while (TRUE) {
        for (j in 1:maxTry) {
          # Correlated random walk
          newAngles <- rwrpnorm(nProposal, angles[step-1], CRWCorrelation)
          #newDistances <- rweibull(nProposal, shape=2, scale=stepSpeedScale * stepIntervalHours * distanceScale)
          newDistances <- randomizeDistance(nProposal)
          
          if (length(homeRangeRadius) != 0 & !is.na(homeRangeRadius)) {
            if (step > 2 & euclidean(coords[1,,drop=F], coords[step-1,,drop=F]) > homeRangeRadius) {
              # Biased correlated random walk
              xy <- matrix(coords[1,,drop=F], ncol=2, nrow=nProposal, byrow=T) - matrix(coords[step-1,,drop=F], ncol=2, nrow=nProposal, byrow=T)
              angleBias <- atan2(xy[,2], xy[,1])
              newAngles <- Arg(BCRWCorrelationBiasTradeoff * exp(complex(imaginary=newAngles)) + (1 - BCRWCorrelationBiasTradeoff) * exp(complex(imaginary=angleBias)))
            }
          }
          
          proposedVectors <- getVector(coords[step-1,,drop=F], newDistances, newAngles)
          acceptedVectors <- randomizeVector(locations=proposedVectors)
          
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
          fileName <- file.path(getwd(), "boundary_reflection_failed_points.RData")
          save(coords, proposedVectors, acceptedVectors, step, file=fileName)
          stop("Boundary reflection failed. File saved to ", fileName)
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
      x <- rep(which(bornIndex), nTransform[bornIndex] - 1)
      survivedBornIndex <- c(which(survivedIndex), x)
      
      agents <<- agents[survivedIndex]
      if (nBorn > 0) agents <<- c(agents, newAgentId:(newAgentId + nBorn - 1))
      newAgentId <<- as.integer(newAgentId + nBorn)
      
      return(list(survivedBornIndex=survivedBornIndex))
    },
    
    randomizeBCRWTracks = function(iteration) {
      library(plyr)
      library(maptools)
      library(rgdal)
      
      initialLocations <- initialPopulation$randomize(nAgents)
      habitatTypes <- extract(study$studyArea$habitat, initialLocations)
      if (any(is.na(habitatTypes))) stop("Invalid initial coordinates.")
      
      nProposal <<- if (inherits(habitatWeights, "uninitializedField") | is.null(habitatWeights)) as.integer(1) else as.integer(10)
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
            function(agentIndex, initialLocations, initialAngles, nAgentsCurrent, isFirst, iteration) {
              ## TODO: UNTESTED CODE
              #if ((inherits(habitatWeights, "uninitializedField") | is.null(habitatWeights) != 0) & loadHabitatRasterInMemory) {
              #  study$studyArea$readRasterIntoMemory()                
              #}
              
              message("Iteration = ", iteration, " / ", nIterations, ", agent (", agents[agentIndex], ") = ", agentIndex, " / ", nAgentsCurrent, ", year = ", year,  " / ", years, ", days = ", days, "...")
              track <- randomizeBCRWTrack(initialLocation=initialLocations[agentIndex,,drop=F],
                                         initialAngle=initialAngles[agentIndex],
                                         isFirst=isFirst,
                                         CRWCorrelation=CRWCorrelation,
                                         BCRWCorrelationBiasTradeoff=BCRWCorrelationBiasTradeoff[iteration],
                                         homeRangeRadius=homeRangeRadius[iteration])
              track$agent <- agents[agentIndex]
              
              return(track)
            },
            initialLocations=initialLocations, initialAngles=initialAngles, nAgentsCurrent=nAgentsCurrent, isFirst=isFirst, iteration=iteration, .parallel=TRUE)
          
          if (year < years) {
            rdReturn <- randomizeBirthDeath()
            survivedBornLastStepIndex <- rdReturn$survivedBornIndex * nSteps
            initialLocations <- as.matrix(track[survivedBornLastStepIndex, c("x","y"), drop=F])
            initialAngles <- track[survivedBornLastStepIndex, c("angle")]
            isFirst <- FALSE
            nAgentsCurrent <- length(agents)
          }        
          
          track$year <- year
          date <- as.POSIXct(strptime(paste(2000+track$year, track$day, track$hour, track$minute, track$second), format="%Y %j %H %M %S"))
          month <- as.POSIXlt(date)$mon + 1
          retainMonths <- c(1,2)
          retainIndex <- month %in% retainMonths
          track <- track[retainIndex,]
          tracks <- rbind(tracks, track)
        }
        
        # TODO: if no agents, randomize new immigrating ones
      }
      
      return(tracks)
    },

    simulateSingle = function(iteration, save=TRUE) {
      tracksDF <- randomizeBCRWTracks(iteration=as.integer(iteration))
      date <- as.POSIXct(strptime(paste(2000+tracksDF$year, tracksDF$day, tracksDF$hour, tracksDF$minute, tracksDF$second), format="%Y %j %H %M %S"))
      tracks <- SimulatedTracks$new(study=study, preprocessData=save, xy=tracksDF[,c("x","y")], id=tracksDF$agent, date=date, iteration=as.integer(iteration))
      return(invisible(tracks))
    },
    
    simulate = function(restartIteration=1, iterationVector=1:nIterations, save=FALSE) {
      stopifnot(restartIteration <= nIterations)
      
      simulatedTracks <- SimulatedTracksCollection$new(study=study)
      for (i in restartIteration:nIterations) {
        message("Iteration ", i, " of ", nIterations, "...")
        tracksDF <- randomizeBCRWTracks(iteration=i)
        date <- as.POSIXct(strptime(paste(2000+tracksDF$year, tracksDF$day, tracksDF$hour, tracksDF$minute, tracksDF$second), format="%Y %j %H %M %S"))
        tracks <- SimulatedTracks$new(study=study, preprocessData=save, xy=tracksDF[,c("x","y")], id=tracksDF$agent, date=date, iteration=i)
        simulatedTracks$addTracks(tracks)
      }
      
      return(invisible(simulatedTracks))
    },
    
    hasHabitatWeights = function() {
      return(!inherits(mss$habitatWeights, "uninitializedField"))
    }
  )
)

# 6 minutes step length for distance correction
MovementSimulationScenarioIntensive <- setRefClass(
  Class = "MovementSimulationScenarioIntensive",
  contains = "MovementSimulationScenario",
  methods = list(
    initialize = function(nAgents=as.integer(50), nIterations=as.integer(1), years=as.integer(1), days=as.integer(60), stepIntervalHours=0.1, runParallel=T, ...) {
      callSuper(years=years, nAgents=nAgents, nIterations=nIterations, days=days, stepIntervalHours=stepIntervalHours, stepSpeedScale=0.5, CRWCorrelation=0.8, runParallel=runParallel, ...)
      return(invisible(.self))
    },
    
    newInstance = function(context, response="Intensive") {
      callSuper()
      study <<- SimulationStudy$new(response=response)$newInstance(context=context)
      initialPopulation <<- RandomInitialPopulation$new(studyArea=study$studyArea)
      return(invisible(.self))
    }
  )
)

# Correlated random walk in a homogeneous landscape, random initial locations
MovementSimulationScenarioA <- setRefClass(
  Class = "MovementSimulationScenarioA",
  contains = "MovementSimulationScenario",
  methods = list(
    initialize = function(nAgents=as.integer(200), nIterations=as.integer(50), years=as.integer(20), days=as.integer(365), stepIntervalHours=4, runParallel=T, ...) {
      callSuper(years=years, nAgents=nAgents, nIterations=nIterations, days=days, stepIntervalHours=stepIntervalHours, stepSpeedScale=0.5, CRWCorrelation=0.8, runParallel=runParallel, ...)
      return(invisible(.self))
    },
    
    newInstance = function(context, response="A") {
      callSuper()
      study <<- SimulationStudy$new(response=response)$newInstance(context=context)
      initialPopulation <<- RandomInitialPopulation$new(studyArea=study$studyArea)
      return(invisible(.self))
    }
  )
)

# Biased correlated random walk in a homogenous landscape, random initial locations
MovementSimulationScenarioB <- setRefClass(
  Class = "MovementSimulationScenarioB",
  contains = "MovementSimulationScenario",
  fields = list(
  ),
  methods = list(
    initialize = function(nAgents=as.integer(200), nIterations=as.integer(50), years=as.integer(20), days=as.integer(365), stepIntervalHours=4, runParallel=T, ...) {
      callSuper(years=years, nAgents=nAgents, nIterations=nIterations, days=days, stepIntervalHours=stepIntervalHours, stepSpeedScale=0.5, runParallel=runParallel, ...)

      pAgentsA <- 0.1
      nAgentsA <- pAgentsA * nAgents
      nAgentsB <- nAgents - nAgentsA
      BCRWCorrelationBiasTradeoff <<- c(rep(NA, nAgentsA), rep(0.3, nAgentsB))
      homeRangeRadius <<- c(rep(NA, nAgentsA), rep(10000, nAgentsB))
      
      return(invisible(.self))
    },
    
    newInstance = function(context, response="B") {
      callSuper()
      study <<- SimulationStudy$new(response=response)$newInstance(context=context)
      initialPopulation <<- RandomInitialPopulation$new(studyArea=study$studyArea)
      return(invisible(.self))
    }
  )
)

###

# Correlated random walk in a homogeneous landscape, random initial locations, fixed distances
MovementSimulationScenarioA.FixedDistances <- setRefClass(
  Class = "MovementSimulationScenarioA.FixedDistances",
  contains = "MovementSimulationScenarioA",
  methods = list(
    randomizeDistance = function(n) {
      rep(stepSpeedScale * stepIntervalHours * distanceScale, n)
    },
    
    newInstance = function(context, response="A.FixedDistances") {
      return(invisible(callSuper(context=context, response=response)))
    }
  )
)


if (F) {


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