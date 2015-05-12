library(devtools)
install_github("statguy/Winter-Track-Counts")
install_github("ropengov/gisfin")
install_github("ropengov/fmi")

library(parallel)
library(doMC)
registerDoMC(cores=round(detectCores()))
source("~/git/Winter-Track-Counts/setup/WTC-Boot.R")
library(WTC)

context <- Context$new(resultDataDirectory=wd.data.results, processedDataDirectory=wd.data.processed, rawDataDirectory=wd.data.raw, scratchDirectory=wd.scratch, figuresDirectory=wd.figures)
study <- FinlandWTCStudy$new(context=context)

if (F) {
  study$response <- "canis.lupus"
  
  ###
  
  tracks <- study$loadTracks()
  habitatWeights <- CORINEHabitatWeights$new(study=study)
  habitatSelection <- tracks$getHabitatPreferences(habitatWeightsTemplate=habitatWeights, nSamples=30, save=F)
  habitatWeights <- CORINEHabitatWeights$new(study=study)$setHabitatSelectionWeights(habitatSelection)
  
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
  sampleIntervals$saveSampleIntervals()
  
  ###
  
  human <- HumanPopulationDensityCovariates$new(study=study)
  weather <- WeatherCovariates$new(study=study, apiKey=fmiApiKey)
  intersections <- study$loadIntersections(predictDistances=FALSE)
  intersections$saveIntersections()
  intersections$setCovariatesId(tag="distance")
  intersections$addCovariates(human, weather)
  intersections$saveCovariates()
  
  ###
  
  intersections <- FinlandWTCIntersections$new(study=study)$loadIntersections()
  intersections$setCovariatesId(tag="density")
  #intersections$loadCovariates()
  habitat <- HabitatSmoothCovariates$new(study=study, scales=2^(1:10) * 62.5)
  intersections$addCovariates(habitat)
  intersections$saveCovariates()
  
  ###
  
  sampleIntervals <- ThinnedMovementSampleIntervals$new(study=study)$setCovariatesId()
  sampleIntervals$loadSampleIntervals()
  
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
    predicted_dist <- exp(intercept - alpha * log(predict_dt));
    
    vector[n_pred] pred_dist_7halfmin;
    vector[n_pred] pred_dist_15min;
    vector[n_pred] pred_dist_30min;
    vector[n_pred] pred_dist_1h;
    vector[n_pred] pred_dist_2h;
    vector[n_pred] pred_dist_4h;
    pred_dist_7halfmin <- exp(intercept + pred_fixed_model_matrix * fixed_effect - alpha * log(7.5/60))
    pred_dist_15min <- exp(intercept + pred_fixed_model_matrix * fixed_effect - alpha * log(15/60))
    pred_dist_30min <- exp(intercept + pred_fixed_model_matrix * fixed_effect - alpha * log(30/60))
    pred_dist_1h <- exp(intercept + pred_fixed_model_matrix * fixed_effect)
    pred_dist_2h <- exp(intercept + pred_fixed_model_matrix * fixed_effect - alpha * log(2))
    pred_dist_4h <- exp(intercept + pred_fixed_model_matrix * fixed_effect - alpha * log(4))
    }
    '
    
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
      
      pred_fixed_model_matrix = pred_fixed_model_matrix,
      n_pred = ncol(pred_fixed_model_matrix)
    ))
    
    estimationResult <<- stan(model_code=intercept_only_code, data=data, iter=iterations, chains=chains)
    estimationResult <<- stan(model_code=fixed_effects_code, data=data, iter=iterations, chains=chains)
    estimationResult <<- stan(model_code=mixed_effects_code, data=data, iter=iterations, chains=chains)
    
    index <- stringr::str_subset(names(estimationResult@sim$samples[[1]]), "predicted_dist")
    predicted_dist <- do.call(c, lapply(estimationResult@sim$samples[[1]][index], mean))
    x <- data.frame(predict_dt, predicted_dist)
    ggplot(x, aes(predict_dt, predicted_dist)) + geom_line() + geom_point(data=movements, aes(dtH, distKm))
    
    return(invisible(.self))
  }
  
  
}
else {
  study$preprocess(fmiApiKey=fmiApiKey) 
}
