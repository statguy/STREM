#setOldClass("lmerMod")

MovementSampleIntervalsPredict <- function(estimationResult, predictionData, randomEffectsFormula=NA) {
  library(lme4)
  distanceKm <- predict(object=estimationResult, newdata=predictionData, REform=randomEffectsFormula)
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
      
      tracksDF <- ld(tracks$tracks)
      tracksDF <- data.frame(tracksDF, breakDownDate(tracksDF$date))

      intervals <<- ddply(tracksDF, .(burst, yday, year), function(x) {
        s <- sum(x$dt, na.rm=T) / 3600
        if (s < 23 | s > 25) return(NULL)
        distKm <- sum(x$dist) / 1e3
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
      intervals$thinId <<- tracks$thinId
      
      return(invisible(.self))
    },
    
    getThinnedTracksSampleIntervals = function(tracks) {
      thinnedTracksCollection <- SimulatedTracksCollection$new(study=study)
      thinnedTracksCollection$addTracks(tracks)
      intervalsList <- list()
      intervalsList[[1]] <- tracks$getSampleIntervals()
      
      maxThins <- 1000
      for (i in 2:maxThins) {
        #thinnedTracks <- thinnedTracksCollection$getTracks(1)$thin(by=i)
        thinnedTracks <- thinnedTracksCollection$getTracks(i-1)$thin(by=2)
        if (is.null(thinnedTracks)) break
        thinnedIntervals <- thinnedTracks$getSampleIntervals()
        
        retainIndex <- thinnedIntervals$intervals$intervalH <= 12
        if (sum(retainIndex) == 0) break
        thinnedIntervals$intervals <- thinnedIntervals$intervals[retainIndex,]
        retainBurst <- unique(thinnedIntervals$intervals[retainIndex, "burst"])
        thinnedTracks$tracks <- thinnedTracks$tracks[burst = retainBurst]

        intervalsList[[i]] <- thinnedIntervals         
        thinnedTracksCollection$addTracks(thinnedTracks)
      }
      if (i == maxThins) stop("Something wrong with thinning.")
      
      intervals <<- ldply(intervalsList, function(x) x$intervals)
      
      return(invisible(thinnedTracksCollection))
    },
    
    getDailyDistanceData = function() {
      if (nrow(intervals) == 0)
        stop("Run get getSampleIntervals() first.")
      return(intervals)      
      
      intervalsSP <- intervals
      coordinates(intervalsSP) <- ~ x+y
      proj4string(intervalsSP) <- study$studyArea$proj4string
      saveCovariates(intervalsSP, save=FALSE)
      dailyDistanceData <- merge(intervals, covariates, sort=F)
      
      return(dailyDistanceData)
    },
    
    fit = function(model=log(distanceKm) ~ log(intervalH) + (1|individualId) + (1|thinId)) {
      library(lme4)  
      
      dailyDistanceData <- getDailyDistanceData()
      estimationResult <<- lmer(model, data=dailyDistanceData)    
      #estimationResult <<- list(C=exp(fixef(result)["(Intercept)"])[1], alpha=-fixef(result)["log(intervalH)"][1])
      
      return(invisible(.self))
    },
        
    predict = function(predictionData=data.frame(intervalH=seq(0, 24, by=0.1)), randomEffectsFormula=NA) {
      distanceKm <- MovementSampleIntervalsPredict(estimationResult=estimationResult, predictionData=predictionData, randomEffectsFormula=randomEffectsFormula)
      return(invisible(distanceKm))
    },
    
    predictDistances = function(model=log(distanceKm) ~ log(intervalH) + (1|individualId) + (1|thinId)) {
      tracks <- study$loadTracks()
      getThinnedTracksSampleIntervals(tracks=tracks)
      fit(model)
      
      # TODO: remove distances > 1
      intersections <- study$loadIntersections()
      dailyDistancePredictions <- intersections$covariates
      dailyDistancePredictions$intervalH <- 1 # TODO: find this
      dailyDistancePredictions$thinId <- 1
      predictedDistances <- 1000 * exp(predict(dailyDistancePredictions))#, ~(1|thinId))
      
      observedDistances <- tracks$getDistances()
      observedDistances <- observedDistances[!is.na(observedDistances)]
      message("Predicted distances = ", mean(predictedDistances, na.rm=T), " ± ", sd(predictedDistances, na.rm=T))
      
      o <- data.frame(Variable="Observed", Value=observedDistances)
      p <- data.frame(Variable="Predicted", Value=predictedDistances)
      distances <- rbind(o, p)
      
      return(invisible(distances))
    },
    
    plotIntervalDistance = function() {
      library(ggplot2)
      p <- ggplot(intervals, aes(intervalH, distanceKm)) +
        geom_point(aes(color=individualId)) +
        ylab("Distance / day (km)") + xlab("Sampling interval (h)") + theme_bw(18)      
      
      if (!inherits(estimationResult, "uninitializedField")) {
        #predictionDate <- estimationResult@frame[,names(fixef(estimationResult))[-1]]        
        #n <- names(fixef(estimationResult))
        
        # TODO
        #predictionData <- data.frame(intervalH=seq(0, 24, by=0.1))
        #predictionData$distanceKm <- predict(predictionData=predictionData)
        #p <- p + geom_line(data=predictionData, aes(intervalH, distanceKm), color="red")
      }
      
      print(p)
    },
    
    plotDistances = function(distances) {
      library(ggplot2)
      p <- ggplot(distances, aes(x=Value, colour=Variable, group=Variable)) + geom_density(fill=NA) +
        xlab("Distance (km)") + ylab("Density")
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
    
    getDailyDistanceData = function() {
      if (nrow(intervals) == 0)
        stop("Run get getSampleIntervals() first.")     
      if (nrow(covariates) == 0) {
        p <- intervals
        coordinates(p) <- ~ x+y
        proj4string(p) <- study$studyArea$proj4string
        saveCovariates(p, save=FALSE)
      }
      
      intervalsSP <- intervals
      coordinates(intervalsSP) <- ~ x+y
      proj4string(intervalsSP) <- study$studyArea$proj4string
      saveCovariates(intervalsSP, save=FALSE)
      dailyDistanceData <- merge(intervals, covariates, sort=F)
      
      return(dailyDistanceData)
    },
    
    #fit = function(model=log(distanceKm) ~ populationDensity + rrday + snow + tday + log(intervalH) + (1|individualId) + (1|thinId)) {
    #  return(invisible(callSuper(model=model)))
    #},
    
    predictDistances = function(model=log(distanceKm) ~ sqrt(populationDensity) + rrday + snow + tday + log(intervalH) + (1|individualId) + (1|thinId)) {
      loadCovariates()
      callSuper(model=model)
    }      
  )
)
