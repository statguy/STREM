library(devtools)
install_github("statguy/Winter-Track-Counts")
install_github("ropengov/gisfin")
install_github("ropengov/fmi")

library(parallel)
library(doMC)
registerDoMC(cores=round(detectCores()))
source("~/git/Winter-Track-Counts/setup/WTC-Boot.R")
library(WTC)


preprocessIntersections <- function(response) {
  context <- Context$new(resultDataDirectory=wd.data.results, processedDataDirectory=wd.data.processed, rawDataDirectory=wd.data.raw, scratchDirectory=wd.scratch, figuresDirectory=wd.figures)
  study <- FinlandWTCStudy$new(context=context)
  study$response <- response

  intersections <- study$loadIntersections(predictDistances=FALSE)
  intersections$saveIntersections()
}

preprocessHabitatPreferences <- function(response) {
  context <- Context$new(resultDataDirectory=wd.data.results, processedDataDirectory=wd.data.processed, rawDataDirectory=wd.data.raw, scratchDirectory=wd.scratch, figuresDirectory=wd.figures)
  study <- FinlandWTCStudy$new(context=context)
  study$response <- response

  tracks <- study$loadTracks()
  habitatWeights <- CORINEHabitatWeights$new(study=study)
  habitatSelection <- tracks$getHabitatPreferences(habitatWeightsTemplate=habitatWeights, nSamples=30, save=T)
}

preprocessSampleIntervals <- function(response) {
  context <- Context$new(resultDataDirectory=wd.data.results, processedDataDirectory=wd.data.processed, rawDataDirectory=wd.data.raw, scratchDirectory=wd.scratch, figuresDirectory=wd.figures)
  study <- FinlandWTCStudy$new(context=context)
  study$response <- response

  human <- HumanPopulationDensityCovariates$new(study=study)
  weather <- WeatherCovariates$new(study=study, apiKey=fmiApiKey)
  sampleIntervals$setCovariatesId()$addCovariates(human, weather)
  sampleIntervals$saveCovariates()
  
  human <- HumanPopulationDensityCovariates$new(study=study)
  weather <- WeatherCovariates$new(study=study, apiKey=fmiApiKey)
  intersections <- study$loadIntersections(predictDistances=FALSE)
  intersections$setCovariatesId(tag="distance")
  intersections$addCovariates(human, weather)
  intersections$saveCovariates()
}

estimateMovementDistances <- function(response) {
  context <- Context$new(resultDataDirectory=wd.data.results, processedDataDirectory=wd.data.processed, rawDataDirectory=wd.data.raw, scratchDirectory=wd.scratch, figuresDirectory=wd.figures)
  study <- FinlandWTCStudy$new(context=context)
  study$response <- response
  
  covariatesFormula <- ~ rrday + snow + tday -1
  iterations <- 1000
  chains <- 1
  
  library(rstan)
  set_cppo("fast")
  
  mixed_effects_code <- '
    data {
    int<lower=1> n_obs;
    vector[n_obs] distKm;
    vector[n_obs] dtH;
    int<lower=1> n_individuals;
    matrix[n_obs,n_individuals] individual_model_matrix;
    int<lower=1> n_samples;
    matrix[n_obs,n_samples] sample_model_matrix;
    int<lower=1> n_fixed;
    matrix[n_obs,n_fixed] fixed_model_matrix;
    int<lower=1> n_predicts;
    vector[n_predicts] predict_dt;
    
    int<lower=1> n_pred_observations;
    int<lower=1> n_pred_covariates;
    matrix[n_pred_observations,n_pred_covariates] pred_fixed_model_matrix;
    }
    parameters {
    real<lower=0> alpha;
    real intercept;
    real<lower=0> error_variance;
    vector[n_individuals] individual_effect;
    real<lower=0> individual_variance;
    vector[n_samples] sample_effect;
    real<lower=0> sample_variance;
    vector[n_fixed] fixed_effect;
    }
    model {
    vector[n_obs] linear_predictor;
    individual_effect ~ normal(0, individual_variance);
    sample_effect ~ normal(0, sample_variance);
    linear_predictor <- intercept + fixed_model_matrix * fixed_effect + individual_model_matrix * individual_effect + sample_model_matrix * sample_effect - alpha * log(dtH);
    log(distKm) ~ normal(linear_predictor, error_variance);
    }
    generated quantities {
    vector[n_predicts] predicted_dist;
    
    vector[n_pred_observations] pred_dist_7halfmin;
    vector[n_pred_observations] pred_dist_15min;
    vector[n_pred_observations] pred_dist_30min;
    vector[n_pred_observations] pred_dist_1h;
    vector[n_pred_observations] pred_dist_2h;
    vector[n_pred_observations] pred_dist_4h;
    
    predicted_dist <- exp(intercept - alpha * log(predict_dt));
    
    pred_dist_7halfmin <- exp(intercept + pred_fixed_model_matrix * fixed_effect - alpha * log(0.125));
    pred_dist_15min <- exp(intercept + pred_fixed_model_matrix * fixed_effect - alpha * log(0.25));
    pred_dist_30min <- exp(intercept + pred_fixed_model_matrix * fixed_effect - alpha * log(0.50));
    pred_dist_1h <- exp(intercept + pred_fixed_model_matrix * fixed_effect);
    pred_dist_2h <- exp(intercept + pred_fixed_model_matrix * fixed_effect - alpha * log(2));
    pred_dist_4h <- exp(intercept + pred_fixed_model_matrix * fixed_effect - alpha * log(4));
    }
    '
    
  sampleIntervals <- ThinnedMovementSampleIntervals$new(study=study)$setCovariatesId()
  sampleIntervals$loadCovariates()
  movements <- sampleIntervals$aggregate()
  movements$distKm <- movements$dist / 1000
  movements$dtH <- movements$dt / 3600
  movements$burst <- as.factor(paste0(movements$burst, movements$yday))
  movements$individual <- as.factor(as.integer(movements$id))
  
  fixed_model_matrix <- model.matrix(covariatesFormula, movements)
  predict_dt <- seq(0, 12, by=0.1)
  intersections <- FinlandWTCIntersections$new(study=study)$setCovariatesId(tag="distance")$loadCovariates()
  pred_fixed_model_matrix <- model.matrix(covariatesFormula, intersections$intersections)
  
  data <- with(movements, list(
    n_obs = length(distKm),
    distKm = distKm,
    dtH = dtH,
    n_individuals = length(unique(individual)),
    individual_model_matrix = model.matrix(~individual-1, movements),
    n_samples = length(unique(burst)),
    sample_model_matrix = model.matrix(~burst-1, movements),
    predict_dt = predict_dt,
    n_predicts = length(predict_dt),
    n_fixed = ncol(fixed_model_matrix),
    fixed_model_matrix = fixed_model_matrix,
    
    n_pred_observations = nrow(pred_fixed_model_matrix),
    n_pred_covariates = ncol(pred_fixed_model_matrix),
    pred_fixed_model_matrix = pred_fixed_model_matrix      
  ))
  
  estimationResult <- stan(model_code=mixed_effects_code, data=data, iter=iterations, chains=chains)
  distanceEstimatesFileName <- context$getFileName(dir=study$context$scratchDirectory, name="DistanceEstimates", response=study$response, region=study$studyArea$region)
  save(estimationResult, file=distanceEstimatesFileName)
  
  summary(rstan::extract(estimationResult, "fixed_effect")[[1]])
  distance30min <- colMeans(rstan::extract(estimationResult, "pred_dist_30min")[[1]])
  distance1h <- colMeans(rstan::extract(estimationResult, "pred_dist_1h")[[1]])
  distance2h <- colMeans(rstan::extract(estimationResult, "pred_dist_2h")[[1]])
  distance4h <- colMeans(rstan::extract(estimationResult, "pred_dist_4h")[[1]])
  
  summary(distance30min[distance30min < 50])
  summary(distance1h[distance1h < 50])
  summary(distance2h[distance2h < 50])
  summary(distance4h[distance4h < 50])
  summary(movements$human)
  #summary(intersections$intersections$human)
  
  intersections <- FinlandWTCIntersections$new(study=study)$setCovariatesId(tag="density")$loadCovariates()
  intersections$setDistance(distance1h)$saveCovariates()
}

preprocessDensityCovariates <- function(response) {
  context <- Context$new(resultDataDirectory=wd.data.results, processedDataDirectory=wd.data.processed, rawDataDirectory=wd.data.raw, scratchDirectory=wd.scratch, figuresDirectory=wd.figures)
  study <- FinlandWTCStudy$new(context=context)
  study$response <- response
  
  intersections <- FinlandWTCIntersections$new(study=study)$loadIntersections()
  intersections$setCovariatesId(tag="density")
  if (intersections$existCovariates()) intersections$loadCovariates()
  habitat <- HabitatSmoothCovariates$new(study=study, scales=2^(0:10) * 62.5)
  elevation <- ElevationSmoothCovariates$new(study=study, scales=2^(0:6) * 1000)
  human <- HumanDensitySmoothCovariates$new(study=study, scales=2^(0:6) * 1000)
  intersections$addCovariates(habitat, elevation, human)
  intersections$saveCovariates()
}

response <- "canis.lupus"
preprocessIntersections(response)
preprocessHabitatPreferences(response)
preprocessSampleIntervals(response)
estimateMovementDistances(response)
preprocessDensityCovariates(response)




if (F) {
  study$response <- "canis.lupus"
  
  intersections <- study$loadIntersections(predictDistances=FALSE)
  intersections$saveIntersections()
  
  ###
  
  tracks <- study$loadTracks()
  habitatWeights <- CORINEHabitatWeights$new(study=study)
  habitatSelection <- tracks$getHabitatPreferences(habitatWeightsTemplate=habitatWeights, nSamples=30, save=T)
  #habitatWeights <- CORINEHabitatWeights$new(study=study)$setHabitatSelectionWeights(habitatSelection)
  
  ###
  
  tracks <- study$loadTracks()
  sampleIntervals <- ThinnedMovementSampleIntervals$new(study=study)
  sampleIntervals$findSampleIntervals(tracks=tracks$tracks)
  #si <- sampleIntervals$aggregate()
  #ggplot(si, aes(dt/3600, dist)) + geom_point()

  human <- HumanPopulationDensityCovariates$new(study=study)
  weather <- WeatherCovariates$new(study=study, apiKey=fmiApiKey)
  #sampleIntervals$intervals <- sampleIntervals$intervals[1:10,]
  sampleIntervals$setCovariatesId()$addCovariates(human, weather)
  sampleIntervals$saveCovariates()
  
  ###
  
  human <- HumanPopulationDensityCovariates$new(study=study)
  weather <- WeatherCovariates$new(study=study, apiKey=fmiApiKey)
  intersections <- study$loadIntersections(predictDistances=FALSE)
  intersections$setCovariatesId(tag="distance")
  intersections$addCovariates(human, weather)
  intersections$saveCovariates()
  
  ###
  
  intersections <- FinlandWTCIntersections$new(study=study)$loadIntersections()
  intersections$setCovariatesId(tag="density")
  habitat <- HabitatSmoothCovariates$new(study=study, scales=64000)
  intersections$addCovariates(habitat)
  intersections$saveCovariates()
  habitat <- HabitatSmoothCovariates$new(study=study, scales=32000)
  intersections$addCovariates(habitat)
  intersections$saveCovariates()

  habitat <- HabitatSmoothCovariates$new(study=study, scales=125)
  intersections$addCovariates(habitat)
  
  
  habitat <- HabitatSmoothCovariates$new(study=study, scales=2^(0:10) * 62.5)
  intersections$addCovariates(habitat)
  intersections$saveCovariates()
  
  
  ## DEBUG
  intersections <- FinlandWTCIntersections$new(study=study)$loadIntersections()
  intersections$setCovariatesId(tag="density")
  intersections$intersections <- intersections$intersections[990:1001,]
  habitat <- HabitatSmoothCovariates$new(study=study, scales=c(125))
  intersections$addCovariates(habitat)
  intersections$intersections@data
  
  which(coordinates(intersections$intersections) == c(3525000,7038000))
  
  # TODO: NaN values???
  intersections$intersections <- intersections$intersections[16975,]
  intersections$intersections@coords <- matrix(c(3693000,6900000),ncol=2)
  habitat <- HabitatSmoothCovariates$new(study=study, scales=10000)
  intersections$addCovariates(habitat)
  
  intersections$intersections <- intersections$intersections[16975,]
  intersections$intersections@coords <- matrix(c(3693000,6905950),ncol=2)
  habitat <- HabitatSmoothCovariates$new(study=study, scales=100)
  intersections$addCovariates(habitat)
  
  ###
  
  
  fit <- function(covariatesFormula=~ human + rrday + snow + tday -1, iterations=1000, chains=1) {
    library(rstan)
    set_cppo("fast")
    
    intercept_only_code <- '
    data {
    int<lower=1> n_obs;
    vector[n_obs] distKm;
    vector[n_obs] dtH;
    }
    parameters {
    real<lower=0> alpha;
    real intercept;
    real<lower=0> error_variance;
    }
    model {
    vector[n_obs] linear_predictor;
    linear_predictor <- intercept - alpha * log(dtH);
    log(distKm) ~ normal(linear_predictor, error_variance);
    }
    '
    
    fixed_effects_code <- '
    data {
    int<lower=1> n_obs;
    vector[n_obs] distKm;
    vector[n_obs] dtH;
    int<lower=1> n_fixed;
    matrix[n_obs,n_fixed] fixed_model_matrix;
    int<lower=1> n_predicts;
    vector[n_predicts] predict_dt;
    }
    parameters {
    real alpha;
    real intercept;
    real<lower=0> error_variance;
    vector[n_fixed] fixed_effect;
    }
    model {
    vector[n_obs] linear_predictor;
    linear_predictor <- intercept + fixed_model_matrix * fixed_effect - alpha * log(dtH);
    log(distKm) ~ normal(linear_predictor, error_variance);
    }
    generated quantities {
    vector[n_predicts] predicted_dist;
    predicted_dist <- exp(intercept - alpha * log(predict_dt));
    }    
    '
    
    mixed_effects_code <- '
    data {
    int<lower=1> n_obs;
    vector[n_obs] distKm;
    vector[n_obs] dtH;
    int<lower=1> n_individuals;
    matrix[n_obs,n_individuals] individual_model_matrix;
    int<lower=1> n_samples;
    matrix[n_obs,n_samples] sample_model_matrix;
    int<lower=1> n_fixed;
    matrix[n_obs,n_fixed] fixed_model_matrix;
    int<lower=1> n_predicts;
    vector[n_predicts] predict_dt;
    
    int<lower=1> n_pred_observations;
    int<lower=1> n_pred_covariates;
    matrix[n_pred_observations,n_pred_covariates] pred_fixed_model_matrix;
    }
    parameters {
    real<lower=0> alpha;
    real intercept;
    real<lower=0> error_variance;
    vector[n_individuals] individual_effect;
    real<lower=0> individual_variance;
    vector[n_samples] sample_effect;
    real<lower=0> sample_variance;
    vector[n_fixed] fixed_effect;
    }
    model {
    vector[n_obs] linear_predictor;
    individual_effect ~ normal(0, individual_variance);
    sample_effect ~ normal(0, sample_variance);
    linear_predictor <- intercept + fixed_model_matrix * fixed_effect + individual_model_matrix * individual_effect + sample_model_matrix * sample_effect - alpha * log(dtH);
    log(distKm) ~ normal(linear_predictor, error_variance);
    }
    generated quantities {
    vector[n_predicts] predicted_dist;
    
    vector[n_pred_observations] pred_dist_7halfmin;
    vector[n_pred_observations] pred_dist_15min;
    vector[n_pred_observations] pred_dist_30min;
    vector[n_pred_observations] pred_dist_1h;
    vector[n_pred_observations] pred_dist_2h;
    vector[n_pred_observations] pred_dist_4h;

    predicted_dist <- exp(intercept - alpha * log(predict_dt));

    pred_dist_7halfmin <- exp(intercept + pred_fixed_model_matrix * fixed_effect - alpha * log(0.125));
    pred_dist_15min <- exp(intercept + pred_fixed_model_matrix * fixed_effect - alpha * log(0.25));
    pred_dist_30min <- exp(intercept + pred_fixed_model_matrix * fixed_effect - alpha * log(0.50));
    pred_dist_1h <- exp(intercept + pred_fixed_model_matrix * fixed_effect);
    pred_dist_2h <- exp(intercept + pred_fixed_model_matrix * fixed_effect - alpha * log(2));
    pred_dist_4h <- exp(intercept + pred_fixed_model_matrix * fixed_effect - alpha * log(4));
    }
    '
    
    sampleIntervals <- ThinnedMovementSampleIntervals$new(study=study)$setCovariatesId()
    sampleIntervals$loadCovariates()
    movements <- sampleIntervals$aggregate()
    movements$distKm <- movements$dist / 1000
    movements$dtH <- movements$dt / 3600
    movements$burst <- as.factor(paste0(movements$burst, movements$yday))
    movements$individual <- as.factor(as.integer(movements$id))
    
    #subset(movements, burst=="Kira.1.03.142")
    #subset(sampleIntervals$intervals, burst=="Kira.1.03.1" & yday==42)
    
    plot(movements$dtH, movements$distKm)
    ggplot(movements, aes(log(dtH), log(distKm))) + geom_point() + stat_smooth(method="lm")
    
    fixed_model_matrix <- if (missing(covariatesFormula)) matrix(0, nrow=nrow(movements), ncol=1)
    else model.matrix(covariatesFormula, movements)
    #fixed_model_matrix <- model.matrix(covariatesFormula, movements)
    predict_dt <- seq(0, 12, by=0.1)
    
    intersections <- FinlandWTCIntersections$new(study=study)$setCovariatesId(tag="distance")$loadCovariates()
    pred_fixed_model_matrix <- model.matrix(covariatesFormula, intersections$intersections)
    
    data <- with(movements, list(
      n_obs = length(distKm),
      distKm = distKm,
      dtH = dtH,
      n_individuals = length(unique(individual)),
      individual_model_matrix = model.matrix(~individual-1, movements),
      n_samples = length(unique(burst)),
      sample_model_matrix = model.matrix(~burst-1, movements),
      predict_dt = predict_dt,
      n_predicts = length(predict_dt),
      n_fixed = ncol(fixed_model_matrix),
      fixed_model_matrix = fixed_model_matrix,
      
      n_pred_observations = nrow(pred_fixed_model_matrix),
      n_pred_covariates = ncol(pred_fixed_model_matrix),
      pred_fixed_model_matrix = pred_fixed_model_matrix      
    ))
    
    estimationResult <<- stan(model_code=intercept_only_code, data=data, iter=iterations, chains=chains)
    estimationResult <<- stan(model_code=fixed_effects_code, data=data, iter=iterations, chains=chains)
    estimationResult <<- stan(model_code=mixed_effects_code, data=data, iter=iterations, chains=chains)
    
    extract(estimationResult, "fixed_effect")
    
    index <- stringr::str_subset(names(estimationResult@sim$samples[[1]]), "predicted_dist")
    predicted_dist <- do.call(c, lapply(estimationResult@sim$samples[[1]][index], mean))
    x <- data.frame(predict_dt, predicted_dist)
    ggplot(x, aes(predict_dt, predicted_dist)) + geom_line() + geom_point(data=movements, aes(dtH, distKm))
    
    return(invisible(.self))
  }
  
  
}
