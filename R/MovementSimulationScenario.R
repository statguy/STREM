randomizeVector = function(locations, habitat, habitatWeights, boundary) {
  library(sp)
  library(raster)
  
  err <- try({
  if (nrow(locations) < 1 | ncol(locations) != 2) {
    print(locations)
    stop("Invalid locations argument.")
  }  
  
  point <- SpatialPoints(locations[,,drop=F], proj4string=boundary@proj4string)
  
  if (inherits(point, "numeric")) {
    print(locations)
    print(point)
    stop("Invalid locations argument.")
  }
  
  outsideBoundary <- is.na(over(point, boundary))
  if (class(outsideBoundary) == "matrix") outsideBoundary <- outsideBoundary[,1]
  }); if (inherits(err, "try-error")) { message(err); stop("randomizeVector(); err = 1, msg = ", err[1]) }
  
  err <- try({
  if (inherits(habitatWeights, "uninitializedField") | is.null(habitatWeights)) {
    if (all(outsideBoundary == TRUE)) return(NULL)
    return(list(index=1, coords=locations[1,,drop=F]))
  }
  else {
    habitatTypes <- extract(habitat, locations)
    w <- habitatWeights$getWeights(habitatTypes)
    w[outsideBoundary] <- 0
    
    if (all(w==0)) return(NULL)
    k <- sample(1:nrow(locations), 1, prob=w)
    return(list(index=k, coords=locations[k,,drop=F]))
  }
  }); if (inherits(err, "try-error")) { message(err); stop("randomizeVector(); err = 2, msg = ", err[1]) }
}

randomizeBCRWTrack <- function(initialLocation, initialAngle, isFirst, nProposal, habitat, habitatWeights, boundary, CRWCorrelation, BCRWCorrelationBiasTradeoff, homeRangeRadius, days, stepIntervalHours, nSteps, distanceScale, stepSpeedScale) {
  library(CircStats)
  library(sp)
  
  err <- try({
  maxTry <- 10000
  
  coords <- matrix(NA, nrow=nSteps + 1, ncol=2)
  coords[1,] <- initialLocation
  angles <- numeric(nSteps + 1)
  angles[1] <- initialAngle
  
  step <- 2
  }); if (inherits(err, "try-error")) { message(err); stop("randomizeBCRWTrack(); err = 1, msg = ", err[1]) }
  
  while (TRUE) {
    for (j in 1:maxTry) {
      
      err <- try({
      # Correlated random walk
      newAngles <- rwrpnorm(nProposal, angles[step-1], CRWCorrelation)
      newDistances <- rweibull(nProposal, shape=2, scale=stepSpeedScale * stepIntervalHours * distanceScale)
      
      if (length(homeRangeRadius) != 0) {
        if (step > 2 & euclidean(coords[1,,drop=F], coords[step-1,,drop=F]) > homeRangeRadius) {
          # Biased correlated random walk
          xy <- matrix(coords[1,,drop=F], ncol=2, nrow=nProposal, byrow=T) - matrix(coords[step-1,,drop=F], ncol=2, nrow=nProposal, byrow=T)
          angleBias <- atan2(xy[,2], xy[,1])
          newAngles <- Arg(BCRWCorrelationBiasTradeoff * exp(complex(imaginary=newAngles)) + (1 - BCRWCorrelationBiasTradeoff) * exp(complex(imaginary=angleBias)))
        }
      }
      
      proposedVectors <- getVector(coords[step-1,,drop=F], newDistances, newAngles)
      acceptedVectors <- randomizeVector(locations=proposedVectors, habitat=habitat, habitatWeights=habitatWeights, boundary=boundary)
      }); if (inherits(err, "try-error")) { message(err); stop("randomizeBCRWTrack(); err = 2, msg = ", err[1]) }
      
      err <- try({
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
      }); if (inherits(err, "try-error")) { message(err); stop("randomizeBCRWTrack(); err = 3, msg = ", err[1]) }
    }
    
    if (j == maxTry) {
      err <- try({
      #track <- SpatialLines(list(Lines(list(Line(coords[1:step-1,])), ID=1)), proj4string=boundary@proj4string)
      #plot(track, col="blue")
      #plot(study$studyArea$boundary, add=T)
      #points(proposedVectors, col="red", pch="+")
      #points(acceptedVectors$coords, col="green", pch="+")
      fileName <- file.path(getwd(), "boundary_reflection_failed_points.RData")
      save(coords, proposedVectors, acceptedVectors, step, file=fileName)
      stop("Boundary reflection failed. File saved to ", fileName)
      }); if (inherits(err, "try-error")) { message(err); stop("randomizeBCRWTrack(); err = 4, msg = ", err[1]) }
    }
    
    if (step == nSteps + 1) break
    step <- step + 1
  }
  
  err <- try({
  index <- if (isFirst) 1:nSteps else 1:nSteps+1
  stepDays <- rep(1:days, each=24 / stepIntervalHours)
  stepHours <- rep(seq(0, 24-stepIntervalHours, by=stepIntervalHours), days)
  stepMinutes <- (stepHours - trunc(stepHours)) * 60
  stepSeconds <- (stepMinutes - trunc(stepMinutes)) * 60
  return(data.frame(x=coords[index,1], y=coords[index,2], angle=angles[index], day=stepDays, hour=trunc(stepHours), minute=trunc(stepMinutes), second=round(stepSeconds)))
  }); if (inherits(err, "try-error")) { message(err); stop("randomizeBCRWTrack(); err = 5, msg = ", err[1]) }
}

# Combined birth-death process:
# 0 born = individual dies/emigrates
# 1 born = individual survives to the next year
# 2 born = 1 parent survives + 1 offspring/immigration
# etc.
randomizeBirthDeath <- function(param=list(mean=1, sd=1.1), agents, newAgentId) {
  err <- try({
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
  
  agents <- agents[survivedIndex]
  if (nBorn > 0) agents <- c(agents, newAgentId:(newAgentId + nBorn - 1))
  newAgentId <- as.integer(newAgentId + nBorn)
  }); if (inherits(err, "try-error")) { message(err); stop("randomizeBirthDeath(); err = 1, msg = ", err[1]) }
  
  return(list(survivedBornIndex=survivedBornIndex, agents=agents, newAgentId=newAgentId))
}

randomizeBCRWTracks <- function(iteration, nIterations, initialLocations, habitat, habitatWeights, nAgents, boundary, CRWCorrelation, BCRWCorrelationBiasTradeoff, homeRangeRadius, days, years, stepIntervalHours, nSteps, distanceScale, stepSpeedScale) {
  library(plyr)
  library(maptools)
  
  options(error=recover)
  
  err <- try({
  #initialLocations <- initialPopulation$randomize(nAgents)
  habitatTypes <- extract(habitat, initialLocations)
  if (any(is.na(habitatTypes))) stop("Invalid initial coordinates.")
  
  nProposal <- if (inherits(habitatWeights, "uninitializedField") | is.null(habitatWeights)) 1 else 10
  message("Number of proposals = ", nProposal)
  
  agents <<- 1:nAgents
  newAgentId <<- as.integer(nAgents + 1)
  
  tracks <- data.frame()
  initialLocations <- coordinates(initialLocations)
  initialAngles <- runif(nAgents, 0, 2*pi)
  isFirst <- TRUE
  
  nAgentsCurrent <- nAgents
  }); if (inherits(err, "try-error")) { message(err); stop("randomizeBCRWTracks(); err = 1, msg = ", err[1]) }
  
  for (year in 1:years) {
    if (nAgentsCurrent > 0) {
      track <- ldply(1:nAgentsCurrent,
        function(agentIndex, initialLocations, initialAngles, agents, nAgentsCurrent, isFirst, nProposal, iteration, nIterations, habitat, habitatWeights, boundary, CRWCorrelation, BCRWCorrelationBiasTradeoff, homeRangeRadius, days, stepIntervalHours, nSteps, distanceScale, stepSpeedScale) {
          ## TODO: UNTESTED CODE
          #if (length(habitatWeights) != 0 & loadHabitatRasterInMemory) {
          #  study$studyArea$readRasterIntoMemory()                
          #}
          
          err <- try({
          message("Iteration = ", iteration, " / ", nIterations, ", agent (", agents[agentIndex], ") = ", agentIndex, " / ", nAgentsCurrent, ", year = ", year,  " / ", years, ", days = ", days, "...")
          track <- randomizeBCRWTrack(initialLocation=initialLocations[agentIndex,,drop=F], initialAngle=initialAngles[agentIndex], isFirst=isFirst, nProposal=nProposal, habitat=habitat, habitatWeights=habitatWeights, boundary=boundary, CRWCorrelation=CRWCorrelation, BCRWCorrelationBiasTradeoff=BCRWCorrelationBiasTradeoff, homeRangeRadius=homeRangeRadius, days=days, stepIntervalHours=stepIntervalHours, nSteps=nSteps, distanceScale=distanceScale, stepSpeedScale=stepSpeedScale)
          track$agent <- agents[agentIndex]
          }); if (inherits(err, "try-error")) { message(err); stop("randomizeBCRWTracks(); err = 2, msg = ", err[1]) }
          
          return(track)
        },
        initialLocations=initialLocations, initialAngles=initialAngles, agents=agents, nAgentsCurrent=nAgentsCurrent, isFirst=isFirst, nProposal=nProposal, iteration=iteration, nIterations=nIterations, habitat=habitat, habitatWeights=habitatWeights, boundary=boundary, CRWCorrelation=CRWCorrelation, BCRWCorrelationBiasTradeoff=BCRWCorrelationBiasTradeoff, homeRangeRadius=homeRangeRadius, days=days, stepIntervalHours=stepIntervalHours, nSteps=nSteps, distanceScale=distanceScale, stepSpeedScale=stepSpeedScale, .parallel=TRUE)
      
      track$year <- year
    }
    
    if (year < years) {
      err <- try({
      rdReturn <- randomizeBirthDeath(agents=agents, newAgentId=newAgentId)
      survivedBornLastStepIndex <- rdReturn$survivedBornIndex * nSteps
      agents <- rdReturn$agents
      newAgentId <- rdReturn$newAgentId
      initialLocations <- as.matrix(track[survivedBornLastStepIndex, c("x","y"), drop=F])
      initialAngles <- track[survivedBornLastStepIndex, c("angle")]
      isFirst <- FALSE
      nAgentsCurrent <- length(rdReturn$agents)
      }); if (inherits(err, "try-error")) { message(err); stop("randomizeBCRWTracks(); err = 3, msg = ", err[1]) }
    }
    
    err <- try({
    date <- as.POSIXct(strptime(paste(2000+track$year, track$day, track$hour, track$minute, track$second), format="%Y %j %H %M %S"))
    month <- as.POSIXlt(date)$mon + 1
    retainMonths <- c(1,2)
    retainIndex <- month %in% retainMonths
    track <- track[retainIndex,]
    
    tracks <- rbind(tracks, track)
    }); if (inherits(err, "try-error")) { message(err); stop("randomizeBCRWTracks(); err = 4") }
    
  }
  
  return(tracks)
}

saveSimulatedTracks <- function(xy, id, date, iteration, tracksDir, response, region) {
  tracks <- data.frame(xy, id=id, date=date)
  thinId <- 1
  fileName <- file.path(tracksDir, paste(paste("Tracks", response, region, iteration, sep="-"), ".RData", sep=""))
  save(tracks, iteration, thinId, file=fileName)
}


##################################


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
    BCRWCorrelationBiasTradeoff = "numeric",
    homeRangeRadius = "numeric",
    distanceScale = "numeric",
    loadHabitatRasterInMemory = "logical",
    runParallel = "logical",
    
    nSteps = "integer"
  ),
  methods = list(
    initialize = function(...) {
      distanceScale <<- 1e3
      loadHabitatRasterInMemory <<- FALSE
      callSuper(...)
      return(invisible(.self))
    },
    
    newInstance = function() {
      nSteps.tmp <- 24 * days / stepIntervalHours
      message("Number of steps = ", nSteps.tmp, ", steps per day = ", 24 / stepIntervalHours)
      if (nSteps.tmp %% 1 != 0) stop("Number of steps must be integer.")
      nSteps <<- as.integer(nSteps.tmp)
      return(invisible(.self))
    },

    simulate = function(restartIteration=1, save=FALSE, useCluster=FALSE) {
      stopifnot(restartIteration <= nIterations)
      
      if (useCluster) {
        library(CNPCluster)       

        localEnv <- new.env()
        localEnv$tracksDir <- study$context$scratchDirectory
        localEnv$response <- study$response
        localEnv$region <- study$studyArea$region
        localEnv$initialLocations <- initialLocations <- initialPopulation$randomize(nAgents)
        
        cnpClusterStartRemote(hosts=cnpClusterGetHostsUkko(maxNodes=min(nIterations, 50), blacklist=c("ukko057.hpc.cs.helsinki.fi")))
        cnpClusterExportCNPCluster()

        cnpClusterExport(c("saveSimulatedTracks", "randomizeBCRWTracks", "randomizeBCRWTrack", "randomizeBirthDeath", "getVector"))
        cluster <- cnpClusterGetRemoteCluster()
        clusterExport(cl=cluster, varlist=c("tracksDir", "response", "region", "initialLocations"), envir=localEnv)
        clusterExport(cl=cluster, varlist=c(
            "nIterations", "nAgents", "habitatWeigts", "CRWCorrelation", "BCRWCorrelationBiasTradeoff",
            "homeRangeRadius", "days", "years", "stepIntervalHours", "nSteps", "distanceScale", "stepSpeedScale"),
          envir=as.environment(.self))
        
        cnpClusterListApplyGeneric(1:nIterations, function(i) {
          message("Iteration ", i, " of ", nIterations, "...")
          tracksDF <- randomizeBCRWTracks(iteration=i, nIterations=nIterations, initialLocations=initialLocations,
                                          habitat=study$studyArea$habitat, habitatWeights=habitatWeights, nAgents=nAgents,
                                          boundary=study$studyArea$boundary, CRWCorrelation=CRWCorrelation,
                                          BCRWCorrelationBiasTradeoff=BCRWCorrelationBiasTradeoff, homeRangeRadius=homeRangeRadius,
                                          days=days, years=years, stepIntervalHours=stepIntervalHours, nSteps=nSteps,
                                          distanceScale=distanceScale, stepSpeedScale=stepSpeedScale)
          date <- as.POSIXct(strptime(paste(2000+tracksDF$year, tracksDF$day, tracksDF$hour, tracksDF$minute, tracksDF$second), format="%Y %j %H %M %S"))
          saveSimulatedTracks(xy=tracksDF[,c("x","y")], id=tracksDF$agent, date=date, iteration=i, tracksDir=tracksDir, response=response, region=region)
        })
          
        cnpClusterStopRemote()
      }
      else {
        #simulatedTracks <- SimulatedTracksCollection$new(study=study)
        for (i in restartIteration:nIterations) {
          message("Iteration ", i, " of ", nIterations, "...")
          initialLocations <- initialPopulation$randomize(nAgents)
          tracksDF <- randomizeBCRWTracks(iteration=i, nIterations=nIterations, initialLocations=initialLocations, habitat=study$studyArea$habitat, habitatWeights=habitatWeights, nAgents=nAgents, boundary=study$studyArea$boundary, CRWCorrelation=CRWCorrelation, BCRWCorrelationBiasTradeoff=BCRWCorrelationBiasTradeoff, homeRangeRadius=homeRangeRadius, days=days, years=years, stepIntervalHours=stepIntervalHours, nSteps=nSteps, distanceScale=distanceScale, stepSpeedScale=stepSpeedScale)
          date <- as.POSIXct(strptime(paste(2000+tracksDF$year, tracksDF$day, tracksDF$hour, tracksDF$minute, tracksDF$second), format="%Y %j %H %M %S"))
          tracks <- SimulatedTracks$new(study=study, preprocessData=save, xy=tracksDF[,c("x","y")], id=tracksDF$agent, date=date, iteration=i)
          #if (returnTracks) simulatedTracks$addTracks(tracks)
        }
        #return(invisible(simulatedTracks))
      }
      
      return(invisible(.self))
    }        
  )
)

# 6 minutes step length for distance fix
MovementSimulationScenarioIntensive <- setRefClass(
  Class = "MovementSimulationScenarioIntensive",
  contains = "MovementSimulationScenario",
  methods = list(
    initialize = function(nAgents=as.integer(50), nIterations=as.integer(6*2), years=as.integer(1), days=as.integer(60), stepIntervalHours=1, runParallel=T, ...) {
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

if (F) {

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