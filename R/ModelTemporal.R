SmoothModelTemporal <- setRefClass(
  Class = "SmoothModelTemporal",
  contains = "Model",
  fields = list(
  ),
  methods = list(
    saveEstimates = function(fileName=getEstimatesFileName()) {
      message("Saving result to ", fileName, "...")
      save(locations, data, model, coordsScale, years, result, offsetScale, family, file=fileName)
    },
    
    estimate = function(save=FALSE, fileName=getEstimatesFileName(), verbose=TRUE) {
      library(INLA)
      
      message("Estimating population density...")
      
      result <<- inla(model,
                      family=family,
                      data=data,
                      E=getObservedOffset(),
                      verbose=verbose,
                      control.predictor=list(compute=TRUE),
                      control.compute=list(cpo=FALSE, dic=TRUE, config=TRUE))
      
      if (is.null(result$ok) | result$ok == FALSE) {
        warning("INLA failed to run.")
      }
      else {
        if (save) saveEstimates(fileName=fileName)
      }
    },
    
    collectEstimates = function(observationWeights=1, predictionWeights=1) {
      library(INLA)
      
      index <- 1:nrow(data)
      data$eta <<- result$summary.linear.predictor$mean[index] + log(observationWeights)
      data$fittedMean <<- result$summary.fitted.values$mean[index] * observationWeights
      data$fittedSD <<- result$summary.fitted.values$sd[index] * observationWeights
      observedOffset <- getObservedOffset()
      
      message("Fitted values sums all years:")
      message("observed = ", sum(data$intersections))
      message("estimated = ", sum(data$fittedMean * observedOffset))
      
      #mean(estimates$data$intersections)
      #estimates$data$fittedMean * estimates$getObservedOffset()
      
      return(invisible(.self))
    },
    
    collectHyperparameters = function() {
      library(INLA)
      
      message("Processing hyperparameters...")
      
      x <- data.frame()
      if (any(rownames(result$summary.hyperpar)=="Precision for year")) # RW1, RW2, AR1, ARp, seasonal
        x <- rbind(x, random_effect_precision=result$summary.hyperpar["Precision for year",])
      if (any(rownames(result$summary.hyperpar)=="Rho for year")) # AR1
        x <- rbind(x, rho=result$summary.hyperpar["Rho for year",])
      if (any(rownames(result$summary.hyperpar)=="PACF1 for year")) # ARp
        x <- rbind(x, rho=result$summary.hyperpar["PACF1 for year",])
      if (any(rownames(result$summary.hyperpar)=="PACF2 for year")) # ARp
        x <- rbind(x, rho=result$summary.hyperpar["PACF2 for year",])
      if (any(rownames(result$summary.hyperpar)=="PACF3 for year")) # ARp
        x <- rbind(x, rho=result$summary.hyperpar["PACF3 for year",])
      if (any(rownames(result$summary.hyperpar)=="PACF4 for year")) # ARp
        x <- rbind(x, rho=result$summary.hyperpar["PACF4 for year",])
      
      return(x)
    }
  )
)

SmoothModelMeanTemporal <- setRefClass(
  Class = "SmoothModelMeanTemporal",
  contains = c("AggregatedModel","SmoothModelTemporal"),
  fields = list(
  ),
  methods = list(
    setModelName = function(family, randomEffect) {
      modelName <<- paste("SmoothModelMean", family, randomEffect, sep="-")
      return(invisible(.self))
    },
    
    setupPrecisionPrior = function(priorParams) {
      precPrior <- list(param=c(priorParams$shape, priorParams$rate), initial=priorParams$initial)
      z <- qgamma(c(0.025, 0.5, 0.975), priorParams$shape, priorParams$rate, log.p=TRUE)
      message("Precision prior 0.025 0.5 0.975 quantiles = ", signif(z[1],2), " ", signif(z[2],2), " ", signif(z[3],2))
      return(precPrior)
    },
    
    setupTemporalPrior = function(priorParams) {
      rhoPrior <- list(param=c(priorParams$mean, priorParams$sd), initial=priorParams$initial)
      z <- qnorm(c(0.025, 0.5, 0.975), priorParams$mean, priorParams$sd)
      message("Rho prior 0.025 0.5 0.975 quantiles = ", signif(z[1],2), " ", signif(z[2],2), " ", signif(z[3],2))
      return(rhoPrior)
    },
    
    setup = function(intersections, params) {
      callSuper(intersections, params)      
      .self$aggregate()
      return(invisible(.self)) 
    },
    
    saveEstimates = function(fileName=getEstimatesFileName()) {
      message("Saving result to ", fileName, "...")
      save(locations, data, model, coordsScale, years, result, offsetScale, family, file=fileName)
    },
    
    estimate = function(save=FALSE, fileName=getEstimatesFileName(), verbose=TRUE) {
      library(INLA)
      
      message("Estimating population density...")
      
      result <<- inla(model,
                      family=family,
                      data=data,
                      E=getObservedOffset(),
                      verbose=verbose,
                      control.predictor=list(compute=TRUE),
                      control.compute=list(cpo=FALSE, dic=TRUE, config=TRUE))
      
      if (is.null(result$ok) | result$ok == FALSE) {
        warning("INLA failed to run.")
      }
      else {
        if (save) saveEstimates(fileName=fileName)
      }
    },
    
    collectEstimates = function(observationWeights=1, predictionWeights=1) {
      library(INLA)
      
      index <- 1:nrow(data)
      data$eta <<- result$summary.linear.predictor$mean[index] + log(observationWeights)
      data$fittedMean <<- result$summary.fitted.values$mean[index] * observationWeights
      data$fittedSD <<- result$summary.fitted.values$sd[index] * observationWeights
      observedOffset <- getObservedOffset()
      
      message("Fitted values sums all years:")
      message("observed = ", sum(data$intersections))
      message("estimated = ", sum(data$fittedMean * observedOffset))      
      
      #mean(estimates$data$intersections)
      #estimates$data$fittedMean * estimates$getObservedOffset()
      
      return(invisible(.self))
    },
    
    collectHyperparameters = function() {
      library(INLA)
      
      message("Processing hyperparameters...")
      
      x <- data.frame()
      if (any(rownames(result$summary.hyperpar)=="Precision for year")) # RW1, RW2, AR1, ARp, seasonal
        x <- rbind(x, random_effect_precision=result$summary.hyperpar["Precision for year",])
      if (any(rownames(result$summary.hyperpar)=="Rho for year")) # AR1
        x <- rbind(x, rho=result$summary.hyperpar["Rho for year",])
      if (any(rownames(result$summary.hyperpar)=="PACF1 for year")) # ARp
        x <- rbind(x, rho=result$summary.hyperpar["PACF1 for year",])
      if (any(rownames(result$summary.hyperpar)=="PACF2 for year")) # ARp
        x <- rbind(x, rho=result$summary.hyperpar["PACF2 for year",])
      if (any(rownames(result$summary.hyperpar)=="PACF3 for year")) # ARp
        x <- rbind(x, rho=result$summary.hyperpar["PACF3 for year",])
      if (any(rownames(result$summary.hyperpar)=="PACF4 for year")) # ARp
        x <- rbind(x, rho=result$summary.hyperpar["PACF4 for year",])
      
      return(x)
    }
  )
)

SimulatedSmoothModelMeanTemporal <- setRefClass(
  Class = "SimulatedSmoothModelMeanTemporal",
  contains = "SmoothModelMeanTemporal",
  fields = list(
    iteration = "integer"  
  ),
  methods = list(
    getEstimatesFileIterations = function() {
      if (inherits(study, "undefinedField") | length(modelName) == 0)
        stop("Provide study and modelName parameters.")
      return(study$context$getIterationIds(dir=study$context$scratchDirectory, name=modelName, response=study$response, region=study$studyArea$region, tag="(\\d+)"))
    },
    
    getEstimatesFileName = function() {
      if (inherits(study, "undefinedField") | length(modelName) == 0 | length(iteration) == 0)
        stop("Provide study, modelName and iteration parameters.")
      return(study$context$getLongFileName(study$context$scratchDirectory, name=modelName, response=study$response, region=study$studyArea$region, tag=iteration))
    }
  )
)

SimulatedSmoothModelTemporal <- setRefClass(
  Class = "SimulatedSmoothModelTemporal",
  contains = "SmoothModelTemporal",
  fields = list(
    iteration = "integer"
  ),
  methods = list(
    getEstimatesFileIterations = function() {
      if (inherits(study, "undefinedField") | length(modelName) == 0)
        stop("Provide study and modelName parameters.")
      return(study$context$getIterationIds(dir=study$context$scratchDirectory, name=modelName, response=study$response, region=study$studyArea$region, tag="(\\d+)"))
    },
    
    getEstimatesFileName = function() {
      if (inherits(study, "undefinedField") | length(modelName) == 0 | length(iteration) == 0)
        stop("Provide study, modelName and iteration parameters.")
      return(study$context$getLongFileName(study$context$scratchDirectory, name=modelName, response=study$response, region=study$studyArea$region, tag=iteration))
    }
  )
)

FinlandSmoothModelTemporal <- setRefClass(
  Class = "FinlandSmoothModelTemporal",
  contains = c("SmoothModelTemporal", "FinlandCovariates"),
  fields = list(
  ),
  methods = list(
    initialize = function(...) {
      callSuper(covariatesName="FinlandSmoothModelCovariates", ...)
      return(invisible(.self))
    },
    
    predictDistances = function(formula=study$getDistanceCovariatesModel(), intervalH=study$getTrackSampleInterval()) {
      distances <- study$predictDistances(formula=formula, data=covariates, intervalH=intervalH)
      return(distances)
    }    
  )
)
