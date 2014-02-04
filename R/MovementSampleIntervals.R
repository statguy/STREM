MovementSampleIntervals <- setRefClass(
  Class = "MovementSampleIntervals",
  fields = list(
    intervals = "data.frame",
    predictions = "data.frame"
  ),
  methods = list(
    initialize = function() {
      callSuper(intervals=data.frame())
    },
    
    # Determines sampling intervals of the recorded movements for each day
    getSampleIntervals = function(tracks) {
      library(plyr)
      
      tracksDF <- ld(tracks$tracks)
      tracksDF$row <- 1:nrow(tracksDF)
      date <- as.POSIXlt(tracksDF$date)
      tracksDF$yday <- date$yday
      tracksDF$year <- date$year
      intervals <<- ddply(tracksDF, .(burst, yday, year), function(x) {
        s <- sum(x$dt, na.rm=T) / 3600
        if (s < 23 | s > 25) return(NULL)
        distKm <- sum(x$dist) / 1e3
        if (is.na(distKm) | distKm > 100) return(NULL)
        intervalMin <- 24 / nrow(x) * 60
        if (intervalMin > 24*60) return(NULL)
        x$intervalSec <- intervalMin * 60
        x$intervalMin <- intervalMin
        x$intervalH <- intervalMin / 60
        x$distanceKm <- distKm
        return(x[1, c("id","date","intervalH","intervalMin","intervalSec","distanceKm","thinid")])
      })
      
      if (nrow(intervals) == 0) warning("Unable to determine sampling intervals.")
    },
    
    getThinnedTracksSampleIntervals = function(tracks) {
      tracksDF <- ld(tracks$tracks)
      date <- as.POSIXlt(tracksDF$date)
      tracksDF$thinid <- as.factor(paste(date$year, date$yday, tracksDF$burst))
      tracks$tracks <- dl(tracksDF)      
      
      thinnedTracksCollection <- SimulatedTracksCollection$new()
      thinnedTracksCollection$add(tracks)
      intervalsList <- list()
      intervalsList[[1]] <- tracks$getSampleIntervals()
      
      maxThins <- 1000
      for (i in 2:maxThins) {
        #thinnedTracks <- thinnedTracksCollection$get(1)$thin(by=i)
        thinnedTracks <- thinnedTracksCollection$get(i-1)$thin(by=2)
        if (is.null(thinnedTracks)) break
        thinnedIntervals <- thinnedTracks$getSampleIntervals()
        
        retainIndex <- thinnedIntervals$intervals$intervalSec <= 60*60*12
        retainIdBurst <- unique(thinnedIntervals$intervals[retainIndex,c("id","burst")])
        if (sum(retainIndex) == 0) break
        thinnedTracks$tracks <- thinnedTracks$tracks[id = retainIdBurst$id][burst = retainIdBurst$burst]
        thinnedIntervals$intervals <- thinnedIntervals$intervals[retainIndex,]
        
        intervalsList[[i]] <- thinnedIntervals         
        thinnedTracksCollection$add(thinnedTracks)
      }
      if (i == maxThins) stop("Something wrong with thinning.")
      
      intervals <<- ldply(intervalsList, function(x) x$intervals)
      
      return(invisible(thinnedTracksCollection))
    },
    
    fit = function() {
      library(lme4)
      result <- lmer(log(distanceKm) ~ log(intervalH) + (1|id) + (1|thinid), data=intervals)
      print(result)
      predictions <<- data.frame(intervalH=seq(0, 24, by=0.1))
      predictions$distanceKm <<- exp(predict(object=result, newdata=predictions, REform=NA))
    },    
    
    plotIntervalDistance = function() {
      library(ggplot2)
      p <- ggplot(intervals, aes(intervalH, distanceKm)) +
        geom_point(aes(color=id)) +
        ylab("Distance / day (km)") + xlab("Sampling interval (h)") + theme_bw(18)      
      if (nrow(predictions) > 0) p <- p + geom_line(data=predictions, aes(intervalH, distanceKm), color="red")
      print(p)
    }
  )
)