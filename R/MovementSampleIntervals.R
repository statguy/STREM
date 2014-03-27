#setOldClass("lmerMod")

MovementSampleIntervalsPredict <- function(estimationResult, predictionData, randomEffectsFormula=NA) {
  if (inherits(estimationResult, "lm")) {
    distanceKm <- predict(estimationResult, newdata=predictionData)
  }
  else {
    library(lme4)
    distanceKm <- predict(object=estimationResult, newdata=predictionData, REform=randomEffectsFormula)
  }
  return(distanceKm)
}

MovementSampleIntervals <- setRefClass(
  Class = "MovementSampleIntervals",
  fields = list(
    study = "Study",
    intervals = "data.frame",
    estimationResult = "ANY"
  ),
  methods = list(
    initialize = function(...) {
      callSuper(...)
      return(invisible(.self))
    },
        
    # Determines sampling intervals of the recorded movements for each day
    getSampleIntervals = function(tracks) {
      library(plyr)
      
      tracksDF <- if (is.data.frame(tracks$tracks)) tracks$tracks else ld(tracks$tracks)
      tracksDF <- data.frame(tracksDF, breakDownDate(tracksDF$date))

      intervals <<- ddply(tracksDF, .(burst, yday, year), function(x) {
        if (is.na(x$dt[nrow(x)])) x$dt[nrow(x)] <- mean(x$dt, na.rm=T)
        
        s <- sum(x$dt, na.rm=T) / 3600
        if (s < 23 | s > 25) return(NULL)
        distKm <- sum(x$dist, na.rm=T) / 1e3
        if (is.na(distKm) | distKm > 100) return(NULL)
        intervalMin <- 24 / nrow(x) * 60
        if (intervalMin > 24*60) return(NULL)
        
        y <- data.frame(id=paste(x$burst[1], x$year[1], x$yday[1]),
                        individualId = x$id[1],
                        date=x$date[1],
                        year=x$year[1],
                        intervalH=intervalMin / 60,
                        intervalMin=intervalMin,
                        intervalSec=intervalMin * 60,
                        distanceKm=distKm,
                        x=mean(x$x),
                        y=mean(x$y))        
        return(y)
      }, .parallel=TRUE)
      
      if (nrow(intervals) == 0) warning("Unable to determine sampling intervals.")
      else intervals$thinId <<- tracks$thinId
      
      return(invisible(.self))
    },
    
    # TODO: rewrite
    getThinnedTracksSampleIntervals = function(tracks, maxThins=200) {
      thinnedTracksCollection <- TracksCollection$new(study=study)
      if (inherits(tracks$tracks, "ltraj")) {
        #tracks$tracks <- ld(tracks$tracks)
        tracks$tracks <- addDtDist(ld(tracks$tracks)) # TODO: remove ltraj objects everywhere
      }
      thinnedTracksCollection$addTracks(tracks)
      intervalsList <- list()
      intervalsList[[1]] <- tracks$getSampleIntervals()
      
      for (i in 2:maxThins) {
        thinnedTracks <- thinnedTracksCollection$getTracks(1)$thin(by=i, thinId=i)
        #thinnedTracks <- thinnedTracksCollection$getTracks(i-1)$thin(by=2, thinId=i)
        #message("nrows before = ", nrow(thinnedTracksCollection$getTracks(1)$tracks), ", nrows after = ", nrow(thinnedTracks$tracks))
        
        if (is.null(thinnedTracks)) break
        
        thinnedIntervals <- thinnedTracks$getSampleIntervals()
        #head(intervalsList[[1]]$intervals) # before
        #head(thinnedIntervals$intervals) # after
        
        if (nrow(thinnedIntervals$intervals) == 0) next        
        
        retainIndex <- thinnedIntervals$intervals$intervalH <= 12
        if (sum(retainIndex) == 0) break

        thinnedIntervals$intervals <- thinnedIntervals$intervals[retainIndex,]
        retainBursts <- unique(thinnedIntervals$intervals[, "burst"])
        thinnedTracks$tracks <- subset(thinnedTracks$tracks, burst %in% retainBursts)
        
        intervalsList[[i]] <- thinnedIntervals  
        thinnedTracksCollection$addTracks(thinnedTracks)
      }
      #if (i == maxThins) stop("Something wrong with thinning.")
      
      intervals <<- ldply(intervalsList, function(x) x$intervals)
      
      return(invisible(thinnedTracksCollection))
    },
    
    getDailyDistanceData = function() {
      if (nrow(intervals) == 0)
        stop("Run get getSampleIntervals() first.")    
      
      #intervalsSP <- intervals
      #coordinates(intervalsSP) <- ~ x+y
      #proj4string(intervalsSP) <- study$studyArea$proj4string
      #saveCovariates(intervalsSP, save=FALSE)
      #dailyDistanceData <- merge(intervals, covariates, sort=F)
      dailyDistanceData <- intervals
      
      return(dailyDistanceData)
    },
    
    fit = function(model=log(distanceKm) ~ log(intervalH) + (1|individualId) + (1|thinId)) {
      library(lme4)  
      dailyDistanceData <- getDailyDistanceData()
      estimationResult <<- lmer(model, data=dailyDistanceData)
      return(invisible(.self))
    },
    
    fitTest = function(model=log(distanceKm) ~ intervalH) {
      dailyDistanceData <- getDailyDistanceData()
      estimationResult <<- glm(model, data=dailyDistanceData, family=gaussian(link="log"))
      return(invisible(.self))
    },
    
    predict = function(predictionData=data.frame(intervalH=seq(0, 12, by=0.1)), randomEffectsFormula=NA) {
      distanceKm <- MovementSampleIntervalsPredict(estimationResult=estimationResult, predictionData=predictionData, randomEffectsFormula=randomEffectsFormula)
      return(invisible(distanceKm))
    },
    
    plotIntervalDistance = function(predictionData=data.frame(intervalH=seq(0, 12, by=0.1)), randomEffectsFormula=NA) {
      library(ggplot2)
      p <- ggplot(intervals, aes(intervalH, distanceKm)) + geom_point() +
        ylab("Distance / day (km)") + xlab("Sampling interval (h)") + theme_bw(18)  
      
      if (!inherits(estimationResult, "uninitializedField")) {
        library(lme4)
        predictionData$distanceKm <- exp(predict(predictionData=predictionData, randomEffectsFormula=randomEffectsFormula))
        p <- p + geom_line(data=predictionData, aes(intervalH, distanceKm), color="red")
      }
      
      print(p)
      return(invisible(p))
    },
    
    plotDistances = function(distances) {
      library(ggplot2)
      p <- ggplot(distances, aes(x=Value, colour=Variable, group=Variable)) + geom_density(fill=NA) +
        xlab("Distance (m)") + ylab("Density")
      print(p)
    }
  )
)

FinlandMovementSampleIntervals <- setRefClass(
  Class = "FinlandMovementSampleIntervals",
  contains = c("MovementSampleIntervals", "FinlandCovariates"),
  methods = list(
    initialize = function(...) {
      callSuper(covariatesName="FinlandMovementSampleIntervalsCovariates", ...)
      return(invisible(.self))
    },
    
    saveIntervalCovariates = function() {
      intervalsSP <- intervals
      coordinates(intervalsSP) <- ~ x+y
      proj4string(intervalsSP) <- study$studyArea$proj4string
      saveCovariates(intervalsSP, save=FALSE)
      return(invisible(.self))
    },
    
    getDailyDistanceData = function() {
      if (nrow(intervals) == 0)
        stop("Run get getSampleIntervals() first.")     
      if (nrow(covariates) == 0) saveIntervalCovariates()
      dailyDistanceData <- merge(intervals, covariates, sort=F)
      
      # To avoid problems with log-transformation
      dailyDistanceData$populationDensity <- dailyDistanceData$populationDensity + 1
      dailyDistanceData$rrday <- dailyDistanceData$rrday + 1
      dailyDistanceData$snow <- dailyDistanceData$snow + 1
      dailyDistanceData$tday <- dailyDistanceData$tday + 100
      
      return(dailyDistanceData)
    },
        
    predictDistances = function(model=log(distanceKm) ~ sqrt(populationDensity) + rrday + snow + tday + log(intervalH) + (1|individualId) + (1|thinId), predictCovariates, truncateDistance=30000) {
      tracks <- study$loadTracks()
      getThinnedTracksSampleIntervals(tracks=tracks)
      fit(model)
      
      # TODO: remove distances > 1
      #intersections <- study$loadIntersections()
      #predictCovariates <- intersections$covariates
      #estimates <- study$loadEstimates()
      #predictCovariates <- estimates$covariates
      
      predictCovariates$intervalH <- 2 # TODO: this needs to be found empirically
      #predictCovariates$thinId <- 1
      predictedDistances <- 1000 * exp(predict(predictCovariates))#, ~(1|thinId))      
      
      message("Estimation results:")
      print(summary(estimationResult))
      
      if (any(predictedDistances > truncateDistance)) {
        message("Invalid predicted distances:")
        predictCovariates$predictedDistances <- predictedDistances
        print(predictCovariates[predictedDistances > truncateDistance,])
        
        # TODO: Extreme predicted distances are truncated for now. Find better prediction model...
        predictedDistances[predictedDistances > truncateDistance] <- truncateDistance
      }
      
      observedDistances <- tracks$getDistances()
      observedDistances <- observedDistances[!is.na(observedDistances)]
      message("Predicted distances = ", mean(predictedDistances, na.rm=T), " ± ", sd(predictedDistances, na.rm=T))
      
      o <- data.frame(Variable="Observed", Value=observedDistances)
      p <- data.frame(Variable="Predicted", Value=predictedDistances)
      distances <- rbind(o, p)
      
      return(invisible(distances))
      #loadCovariates()
      #callSuper(model=model)
    }      
  )
)
