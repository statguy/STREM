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
      library(adehabitatLT)
      library(plyr)
      
      tracksDF <- if (inherits(tracks$tracks, "ltraj")) ld(tracks$tracks) else tracks$tracks
      tracksDF <- data.frame(tracksDF, breakDownDate(tracksDF$date))

      intervals <<- ddply(tracksDF, .(id, yday, year), function(x) {
        if (is.na(x$dt[nrow(x)])) x$dt[nrow(x)] <- mean(x$dt, na.rm=T)
        
        s <- sum(x$dt, na.rm=T) / 3600
        if (s < 23 | s > 25) return(NULL)
        distKm <- sum(x$dist, na.rm=T) / 1e3
        if (is.na(distKm) | distKm > 100) return(NULL)
        intervalMin <- 24 / nrow(x) * 60
        if (intervalMin > 24*60) return(NULL)
        
        y <- data.frame(id=paste(x$id[1], x$year[1], x$yday[1]),
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

FinlandMovementSampleIntervals <- setRefClass(
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
