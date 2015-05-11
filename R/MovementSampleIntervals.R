MovementSampleIntervals <- setRefClass(
  Class = "MovementSampleIntervals",
  fields = list(
    study = "ANY",
    intervals = "data.frame",
    maxDt = "numeric"
  ),
  methods = list(
    findSampleIntervals = function(tracks) {
      stop("Unimplemented method.")
    },
    
    getSampleIntervals = function() intervals,
    
    getSampleIntervalsFileName = function() {
      return(context$getFileName(dir=study$context$processedDataDirectory, name="MovementSampleIntervals", region=study$studyArea$region))
    },
    
    saveSampleIntervals = function(fileName=getSampleIntervalsFileName()) {
      save(intervals, file=fileName)
      return(invisible(.self))
    },
    
    loadSampleIntervals = function(fileName=getSampleIntervalsFileName()) {
      load(fileName, env=as.environment(.self))
      return(invisible(.self))
    }
  )
)

ConstantMovementSampleIntervals <- setRefClass(
  Class = "ConstantMovementSampleIntervals",
  contains = "MovementSampleIntervals",
  fields = list(
    tracks = "ANY"
  ),
  methods = list(
    findSampleIntervals = function(tracks) {
      if (!inherits(tracks, "Tracks"))
        stop("Argument 'tracks' must be of class 'Tracks'.")
      
      .self$tracks <- tracks
      return(invisible(.self))
    },
    
    getSampleIntervals = function() {
      return(tracks$tracks)
    }
  )
)

MaxApproximateConstantMovementSampleIntervals <- setRefClass(
  Class = "MaxApproximateConstantMovementSampleIntervals",
  contains = "MovementSampleIntervals",
  fields = list(
    bandWidth = "numeric"
  ),
  methods = list(
    findSampleIntervals = function(tracks) {
      if (!inherits(tracks, "ltraj"))
        stop("Argument 'tracks' must be of class 'ltraj'.")
      
      if (length(maxDt) == 0) maxDt <<- 24 * 60 * 60
      x <- ld(tracks)
      x <- x[x$dt <= maxDt & !is.na(x$dt),]
      if (length(bandWidth) == 0) bandWidth <<- 10 * 60 # 10 minutes default bandwidth
      d <- density(x$dt, bw=bandWidth, kernel="gaussian")
      retainDt <- d$x[which.max(d$y)] + c(-bandWidth, bandWidth) / 2
      
      message("Sampling interval with max observations: min dt = ", retainDt[1], ", max dt = ", retainDt[2])
      
      intervals <<- x[x$dt >= retainDt[1] & x$dt <= retainDt[2],]
      return(invisible(.self))
    }  
  )
)

.retain24HoursBursts <- function(x) {
  y <- ddply(x, .(burst, yday, thin), function(x) {
    if (all(is.na(x$dt))) return(NULL)
    interpolatedDt <- if (is.na(x$dt[length(x$dt)])) 24 * 3600 / length(x$dt) else 0
    dt <- sum(x$dt, na.rm=T) + interpolatedDt
    if (dt >= 23 * 3600 & dt <= 25 * 3600) return(x)
    return(NULL)
  })
  message(nrow(x)-nrow(y), " of ", nrow(x), " vectors removed due to not summing up to 24 hours.")
  return(y)
}

.retainConstantSamplingIntervalBursts <- function(x) {
  y <- ddply(x, .(burst, yday, thin), function(x) {
    dt <- x$dt[!is.na(x$dt)]
    if (length(dt) <= 1) return(x)
    if (length(x$dt) >= 48 || length(unique(nextn(dt, 2))) == 1) return(x)
    else return(NULL)
  })
  message(nrow(x)-nrow(y), " of ", nrow(x), " vectors removed due to non-constant sampling interval.")
  return(y)
}

.interpolateLastMovement <- function(x) {
  y <- ddply(x, .(burst, yday, thin), function(x) {        
    last <- x[nrow(x),]
    if (is.na(last$dt)) {
      speed <- sum(x$dist, na.rm=T) / sum(x$dt, na.rm=T)
      diffFrom24h <- 24 * 3600 - sum(x$dt, na.rm=T)
      if (diffFrom24h < 0) {
        warning("There are bursts of over 24h with missing last movement.")
        diffFrom24h <- 0
      }
      x[nrow(x),]$dt <- diffFrom24h
      x[nrow(x),]$dist <- diffFrom24h * speed
      #print(x)
    }
    return(x)
  })
  return(y)
}

.retainBurstsWithMovementsLessThan12Hours = function(x, maxDt = 12 * 3600 + .5 * 3600) {
  y <- ddply(x, .(burst, yday, thin), function(x) {
    if (any(x$dt > maxDt)) return(NULL)
    else return(x)
  })
  message(nrow(x)-nrow(y), " of ", nrow(x), " vectors longer than 12 hours removed.")
  return(y)
}

ThinnedMovementSampleIntervals <- setRefClass(
  Class = "ThinnedMovementSampleIntervals",
  contains = "MovementSampleIntervals",
  fields = list(
    maxThinnings = "numeric",
    covariateNames = "character"
  ),
  methods = list(
    getSampleLocations = function() {
      xyt <- subset(intervals, thin == 1)
      sp::coordinates(xyt) <- ~ x+y
      sp::proj4string(xyt) <- study$studyArea$boundary@proj4string
      return(xyt)
    },
        
    findSampleIntervals = function(tracks) {
      library(plyr)
      library(adehabitatLT)
      
      if (!inherits(tracks, "ltraj"))
        stop("Argument 'tracks' must be of class 'ltraj'.")
      
      x <- adehabitatLT::ld(tracks)
      if (length(maxDt) == 0) maxDt <<- 12 * 60 * 60 + 1800 # Remove vectors of over 12 hours (with tolerance)
      #index <- x$dt <= maxDt
      #index[is.na(index)] <- T
      #x <- x[index,]
      x$yday <- as.POSIXlt(x$date)$yday
      x$thin <- 1
      x$observation <- 1:nrow(x)
      
      # Retain movements that cover 24 hours
      y <- .retain24HoursBursts(x)
      
      # Detect if the sampling interval is not constant and remove those bursts, but ignore cases
      # with low intervals (as there is not much low interval data).
      # Since missing movements cause intervals dt to increase by factor of 2 or more, we can find
      # equal or closest number that is power of 2 to dt. If all such numbers are equal, we have no
      # missing movements.
      y <- .retainConstantSamplingIntervalBursts(y)
      
      # 24 hours activity
      #y$hour <- as.POSIXlt(y$date)$hour
      #ggplot(y, aes(hour, dist)) + geom_histogram(stat="identity")
      
      x <- y
      intervals <<- x
      
      if (length(maxThinnings) == 0) maxThinnings <<- 100
      thinFactor <- 2
      for (thinningIndex in 1:maxThinnings) {
        message("Thinning ", thinningIndex, "/", maxThinnings, ", factor = ", thinFactor)
        
        # Thin movements
        y <- plyr::ddply(x, .(burst), function(x, thinFactor) {
          retainIndex <- seq(1, nrow(x), by=thinFactor)
          x <- x[retainIndex,]
          
          x1 <- x[-1,]
          x2 <- x[-nrow(x),]
          x$dist <- c(sqrt((x1$x - x2$x)^2 + (x1$y - x2$y)^2), NA)
          x$dt <- c(unclass(x1$date) - unclass(x2$date), NA)
          
          if (any(x$dt[!is.na(x$dt)] < 0) || any(is.na(x$dt[-length(x$dt)]))) print(x)
          
          index <- x$dt <= maxDt
          index[is.na(index)] <- TRUE
          if (any(index == FALSE)) return(NULL)
          if (nrow(x) <= 1) return(NULL)
          
          x$thin <- thinFactor
          return(x)
        }, thinFactor=thinFactor)
        if (length(y) == 0) break

        intervals <<- rbind(intervals, y)
        thinFactor <- thinFactor + 1
      }
      
      # Retain movements that cover 24 hours (again)
      y <- .retain24HoursBursts(intervals)
      
      # Interpolate the missing last movements by replacing with mean.
      # However, the activity at last hours might be different from daytime and thus the estimates could be biased.
      y <- .interpolateLastMovement(y)
      
      # Remove bursts with movements of over 12 hours (with tolerance)
      intervals <<- .retainBurstsWithMovementsLessThan12Hours(y)

      return(invisible(.self))
    },
    
    associateCovariates = function(...) {
      library(plyr)
      covariates <- cbind(...)
      covariateNames <<- colnames(covariates)
      intervals <<- plyr::ddply(intervals, .(thin), function(x, covariates) {
        return(plyr::join(x, covariates))
      }, covariates=cbind(data.frame(observation=subset(intervals, thin == 1, select="observation")), covariates))
      return(invisible(.self))
    },
    
    addCovariates = function(...) {
      library(plyr)
      covariates <- list(...)
      values <- llply(covariates, function(x) {
        if (!inherits(x, "Covariates"))
          stop("Arguments must be of class 'Covariates'.")
        x$preprocess()
        y <- x$extract(getSampleLocations())
        return(y)
      })      
      do.call(.self$associateCovariates, values)
      return(invisible(.self))
    },
    
    aggregate = function() {
      x <- ddply(intervals, .(burst, yday, thin), function(x) {
        #estDist <- sum(x$dist, na.rm=T) / sum(x$dt, na.rm=T) * 24 * 60 * 60
        estDist <- sum(x$dist, na.rm=T)
        y <- data.frame(id=x$id[1], thin=x$thin[1], year=as.POSIXlt(x$date[1])$year+1900, yday=x$yday[1], dt=mean(x$dt, na.rm=T), dist=estDist,
                   dt.mean=mean(x$dt, na.rm=T), dt.sd=sd(x$dt, na.rm=T))
        meanCovariates <- as.list(colMeans(x[,covariateNames]))
        if (length(meanCovariates) > 0) y <- cbind(y, meanCovariates)
        return(y)
      })
      return(x)      
    }
  )  
)



# DEPRECATED
if (F) {

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

DEPRECATED_MovementSampleIntervals <- setRefClass(
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
      library(adehabitatLT)
      library(dplyr)
      
      tracksDF <- if (inherits(tracks$tracks, "ltraj")) ld(tracks$tracks) else tracks$tracks
      tracksDF <- data.frame(tracksDF, breakDownDate(tracksDF$date))
      
      getInterval <- function(dt, dist) {
        if (is.na(dt[length(dt)])) dt[length(dt)] <- mean(dt, na.rm=T)
        s <- sum(dt, na.rm=T) / 3600
        if (s < 23 | s > 25) return(NA)
        distKm <- sum(dist, na.rm=T) / 1e3
        if (is.na(distKm) | distKm > 100) return(NA)
        intervalMin <- 24 / length(dt) * 60
        if (intervalMin > 24*60) return(NA)
        return(intervalMin)
        #return(list(intervalMin=intervalMin, distKm=distKm))
      }
      
      x <- tracksDF %>% group_by(id, yday, year) %>%
        summarise(id2=paste(id[1], yday[1], year[1]), burst=burst[1], individualId=id[1], date=date[1],
                  intervalMin=getInterval(dt, dist), distanceKm=sum(dist, na.rm=T) / 1e3, x=mean(x), y=mean(y)) %>%
        mutate(intervalH=intervalMin / 60, intervalSec=intervalMin * 60)
      x$id <- x$id2
      x$id2 <- NULL
      x <- x[!is.na(x$intervalMin),]
      intervals <<- as.data.frame(x)
      
      message("Distance distribution:")
      print(summary(x$distanceKm))
      message("Interval distribution:")
      print(summary(x$intervalH))
      message("Speed distribution:")
      print(summary(x$distanceKm / x$intervalH))
      
      if (F) {
      library(plyr)
      intervals <<- ddply(tracksDF, .(id, yday, year), function(x) {
        if (is.na(x$dt[nrow(x)])) x$dt[nrow(x)] <- mean(x$dt, na.rm=T)
        
        s <- sum(x$dt, na.rm=T) / 3600
        if (s < 23 | s > 25) return(NULL)
        distKm <- sum(x$dist, na.rm=T) / 1e3
        if (is.na(distKm) | distKm > 100) return(NULL)
        intervalMin <- 24 / nrow(x) * 60
        if (intervalMin > 24*60) return(NULL)
        
        y <- data.frame(id=paste(x$id[1], x$year[1], x$yday[1]),
                        burst=x$burst[1],
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
      }
      
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
        # TODO: burst definition changed, DOES THIS WORK?
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
    
    fit = function(covariatesFormula, iterations=1000, chains=1) {
      library(rstan)
      set_cppo("fast")
      
      code <- '
      data {
        int<lower=1> n_obs;
        vector[n_obs] distance;
        vector[n_obs] interval;
        int<lower=1> n_individuals;
        matrix[n_obs,n_individuals] individual_model_matrix;
        int<lower=1> n_samples;
        matrix[n_obs,n_samples] sample_model_matrix;
        int<lower=1> n_fixed;
        matrix[n_obs,n_fixed] fixed_model_matrix;
        int<lower=1> n_predicts;
        vector[n_predicts] predict_interval;
      }
      parameters {
        real alpha;
        real intercept;
        real<lower=0> sigma;
        vector[n_individuals] individual_effect;
        real<lower=0> individual_sigma;
        vector[n_samples] sample_effect;
        real<lower=0> sample_sigma;
        vector[n_fixed] fixed_effect;
      }
      model {
        vector[n_obs] mu;
        individual_effect ~ normal(0, individual_sigma);
        sample_effect ~ normal(0, sample_sigma);
        mu <- intercept + fixed_model_matrix * fixed_effect - alpha * log(interval) + individual_model_matrix * individual_effect + sample_model_matrix * sample_effect;
        log(distance) ~ normal(mu, sigma);
      }
      generated quantities {
        vector[n_predicts] predicted_distance;
        predicted_distance <- exp(intercept - alpha * log(predict_interval));
      }
      '
      
      #covariatesFormula <- ~populationDensity + rrday + snow + tday - 1
      movements <- getDailyDistanceData()
      fixed_model_matrix <- if (missing(covariatesFormula)) matrix(0, nrow=nrow(movements), ncol=1) # Hack
      else model.matrix(covariatesFormula, movements)
      
      predict_interval <- seq(0, 12, by=0.1)
      data <- with(movements, list(
        n_obs = length(distanceKm),
        distance = distanceKm,
        interval = intervalH,
        n_individuals = length(unique(individualId)),
        individual_model_matrix = model.matrix(~individualId-1, movements),
        n_samples = length(unique(sampleId)),
        sample_model_matrix = model.matrix(~sampleId-1, movements),
        predict_interval = predict_interval,
        n_predicts = length(predict_interval),
        n_fixed = ncol(fixed_model_matrix),
        fixed_model_matrix = fixed_model_matrix
      ))
      
      estimationResult <<- stan(model_code=code, data=data, iter=iterations, chains=chains)
      #estimationResult@sim$fnames_oi[grep("fixed_effect", estimationResult@sim$fnames_oi)] <<- colnames(fixed_model_matrix)      
      #estimationResult@sim$samples
      
      return(invisible(.self))
    },
    
    predict = function(fixed_model_matrix, intervalH) {
      estimatedValues <- rstan::extract(estimationResult)
      intercept <- mean(estimatedValues$intercept)
      alpha <- mean(estimatedValues$alpha)
      
      if (missing(fixed_model_matrix)) {
        distanceKm <- exp(intercept - alpha * log(intervalH))
      }
      else {
        fixed_effect <- colMeans(estimatedValues$fixed_effect)
        if (ncol(fixed_model_matrix) != length(fixed_effect))
          stop("Model matrix must have the same number of columns as is length of the fixed_effect vector.")
        distanceKm <- exp(intercept + fixed_model_matrix %*% fixed_effect - alpha * log(intervalH))
      }
      
      return(invisible(distanceKm))
    },
    
    getSampleIntervalsFileName = function() {
      return(study$context$getFileName(dir=study$context$resultDataDirectory, "MovementSampleIntervals", response=study$response, region=study$studyArea$region))
    },
    
    saveSampleIntervals = function(fileName=getSampleIntervalsFileName()) {
      save(intervals, estimationResult, file=fileName)
      return(invisible(.self))
    },
    
    loadSampleIntervals = function(fileName=getSampleIntervalsFileName()) {
      load(file=fileName, envir=as.environment(.self))
      return(invisible(.self))
    },
    
    estimatedValuesSummary = function() {
      estimatedValues <- rstan::extract(estimationResult)
      intercept <- summaryStat(estimatedValues$intercept, "intercept")
      alpha <- summaryStat(estimatedValues$alpha, "alpha")
      fixed_effect <- if (any("fixed_effect" == names(estimatedValues)))
        do.call("rbind", apply(estimatedValues$fixed_effect, 2, summaryStat, "fixed_effect"))
      else NULL
      sigma <- summaryStat(estimatedValues$sigma, "sigma")
      individual_sigma <- summaryStat(estimatedValues$individual_sigma, "individual_sigma")
      sample_sigma <- summaryStat(estimatedValues$sample_sigma, "sample_sigma")
      
      print(rbind(intercept, alpha, fixed_effect, sigma, individual_sigma, sample_sigma))
      cat("\n")
      cat(paste("Log likelihood =", mean(estimatedValues$lp__)))
      
      return(invisible(.self))
    },
    
    getDistanceCurve = function(predict_interval=seq(0, 12, by=0.1)) {
      movementsFit <- data.frame(intervalH=predict_interval,
                                 distanceKm=colMeans(rstan::extract(estimationResult)$predicted_distance),
                                 distanceKm25=apply(rstan::extract(estimationResult)$predicted_distance, 2, quantile, .025),
                                 distanceKm975=apply(rstan::extract(estimationResult)$predicted_distance, 2, quantile, .975))
      return(movementsFit)      
    },
    
    plotIntervalDistance = function(plot=FALSE) {
      library(ggplot2)
      
      p <- ggplot(intervals, aes(intervalH, distanceKm)) + geom_point(alpha=.3) +
        ylab("Distance / day (km)") + xlab("Sampling interval (h)")
      
      if (!inherits(estimationResult, "uninitializedField")) {
        movementsFit <- getDistanceCurve()
        p <- p + geom_line(data=movementsFit, aes(intervalH, distanceKm), color="red") +
          geom_smooth(data=movementsFit, aes(ymin=distanceKm25, ymax=distanceKm975), stat="identity") + ylim(c(0, max(intervals$distanceKm)))
      }
      
      if (plot) print(p)
      
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

DEPRECATED_FinlandMovementSampleIntervals <- setRefClass(
  Class = "FinlandMovementSampleIntervals",
  contains = c("MovementSampleIntervals", "FinlandCovariates"),
  methods = list(
    initialize = function(...) {
      callSuper(covariatesName="FinlandMovementSampleIntervalsCovariates", ...)
      return(invisible(.self))
    },
    
    saveTrackCovariates = function(save=FALSE) {
      intervalsSP <- intervals
      coordinates(intervalsSP) <- ~ x+y
      proj4string(intervalsSP) <- study$studyArea$proj4string
      saveCovariates(intervalsSP, save=save)
      if (nrow(covariates) != nrow(intervals))
        stop("Something wrong with covariates.")
      return(invisible(.self))
    },
    
    getDailyDistanceData = function() {
      #library(dplyr)
      
      if (nrow(intervals) == 0)
        stop("Run get getSampleIntervals() first.")     
      if (nrow(covariates) == 0) saveIntervalCovariates()
      #movements <- data.frame(intervals, select(covariates, -c(id,year)))
      movements <- data.frame(intervals, covariates[,!colnames(covariates) %in% c("id","year")])
      movements$sampleId <- with(movements, paste(individualId, date))
      
      # To avoid problems with log-transformation
      #movements$populationDensity <- movements$populationDensity + 1
      #movements$rrday <- movements$rrday + 1
      #movements$snow <- movements$snow + 1
      #movements$tday <- movements$tday + 100
      
      return(movements)
    },
        
    predictDistances = function(model=log(distanceKm) ~ sqrt(populationDensity) + rrday + snow + tday + log(intervalH) + (1|individualId) + (1|thinId), predictCovariates, truncateDistance=30, thin=TRUE) {
      if (thin) {
        tracks <- study$loadTracks()
        getThinnedTracksSampleIntervals(tracks=tracks)
      }
      
      #fit(model) # TODO: fix the problem here
      # quickfix:
      model <- log(distanceKm) ~ sqrt(sqrt(populationDensity)) + sqrt(rrday) + sqrt(snow) + log(intervalH) # problem with tday
      estimationResult <<- glm(model, data=getDailyDistanceData())
      
      predictCovariates$intervalH <- 2 # TODO: this needs to be found empirically
      #predictCovariates$thinId <- 1
      predictedDistances <- exp(predict(predictCovariates))#, ~(1|thinId))
      
      message("Estimation results:")
      print(summary(estimationResult))
      
      if (any(predictedDistances > truncateDistance)) {
        message("Invalid predicted distances:")
        predictCovariates$predictedDistances <- predictedDistances
        print(predictCovariates[predictedDistances > truncateDistance,])
        
        # TODO: Extreme predicted distances are truncated for now. Find better prediction model or add more data, etc...
        predictedDistances[predictedDistances > truncateDistance] <- truncateDistance
      }
      
      #observedDistances <- tracks$getDistances()
      #observedDistances <- observedDistances[!is.na(observedDistances)]
      #message("Predicted distances = ", mean(predictedDistances, na.rm=T), " ± ", sd(predictedDistances, na.rm=T))
      
      #o <- data.frame(Variable="Observed", Value=observedDistances)
      #p <- data.frame(Variable="Predicted", Value=predictedDistances)
      #distances <- rbind(o, p)
      
      return(invisible(predictedDistances))
    }      
  )
)

}
