setOldClass("lmerMod")

MovementSampleIntervals <- setRefClass(
  Class = "MovementSampleIntervals",
  fields = list(
    study = "Study",
    intervals = "data.frame",
    estimationResult = "lmerMod"
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
      #tracksDF$row <- 1:nrow(tracksDF)
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
        if (any("thinId" %in% colnames(x))) y$thinId <- x$thinId[1]
        
        return(y)
      }, .parallel=TRUE)
      
      if (nrow(intervals) == 0) warning("Unable to determine sampling intervals.")
      return(invisible(.self))
    },
    
    getThinnedTracksSampleIntervals = function(tracks) {
      tracksDF <- ld(tracks$tracks)
      date <- as.POSIXlt(tracksDF$date)
      tracksDF$thinId <- as.factor(paste(date$year, date$yday, tracksDF$burst))
      tracks$tracks <- dl(tracksDF)      
      
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
        
        retainIndex <- thinnedIntervals$intervals$intervalSec <= 60*60*12
        retainIdBurst <- unique(thinnedIntervals$intervals[retainIndex,c("id","burst")])
        if (sum(retainIndex) == 0) break
        thinnedTracks$tracks <- thinnedTracks$tracks[id = retainIdBurst$id][burst = retainIdBurst$burst]
        thinnedIntervals$intervals <- thinnedIntervals$intervals[retainIndex,]
        
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
    
    predict = function(predictionData=data.frame(intervalH=seq(0, 24, by=0.1))) {
      library(lme4)
      distanceKm <- exp(predict(object=result, newdata=predictionData, REform=NA))
      return(invisible(distanceKm))
    },
    
    plotIntervalDistance = function() {
      library(ggplot2)
      p <- ggplot(intervals, aes(intervalH, distanceKm)) +
        geom_point(aes(color=id)) +
        ylab("Distance / day (km)") + xlab("Sampling interval (h)") + theme_bw(18)      
      if (nrow(predictions) > 0) p <- p + geom_line(data=predictions, aes(intervalH, distanceKm), color="red")
      print(p)
    },
    
    applyDistanceCorrection = function(surveyRoutes) {
      # TODO
      intervals$distanceKmCorrected <<- estimatedParameters$C * intervals$intervalH ^ -estimatedParameters$alpha
      return(invisible(.self))
    },
    
    getMeanDistance = function() {
      return(mean(intervals$distanceKm))
    }
  )
)

FinlandMovementSampleIntervals <- setRefClass(
  Class = "FinlandMovementSamplingIntervals",
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
    
    fit = function(model=log(distanceKm) ~ populationDensity + rrday + snow + tday + log(intervalH) + (1|individualId) + (1|thinId)) {
      return(invisible(callSuper(model=model)))
    }
  )
)
