# ./parallel_r.py -t 1:3 -n 6 -l 10.0 -b ~/tmp/blacklist.txt -v ~/git/Winter-Track-Counts/inst/wtc/estimate.R notest
# ./parallel_r.py -t 3 -n 2 -l 10.0 -b ~/tmp/blacklist.txt -v ~/git/Winter-Track-Counts/inst/wtc/estimate.R notest

# library(devtools); install_github("statguy/Winter-Track-Counts")

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2) stop("Invalid arguments.")
test <- args[1]
task_id <- args[length(args)]
message("Arguments provided:")
print(args)

library(parallel)
library(doMC)
registerDoMC(cores=detectCores())
library(WTC)
source("~/git/Winter-Track-Counts/setup/WTC-Boot.R")

if (test == "test") {
  # For testing

  context <- Context$new(resultDataDirectory=wd.data.results, processedDataDirectory=wd.data.processed, rawDataDirectory=wd.data.raw, scratchDirectory=wd.scratch, figuresDirectory=wd.figures)
  study <- FinlandWTCStudy$new(context=context, response=response, distanceCovariatesModel=~populationDensity+rrday+snow+tday-1, trackSampleInterval=2)
  
  intersections <- study$loadIntersections()
  model <- FinlandSmoothModelTemporal$new(study=study)
  model$setup(intersections=intersections, params=list(family="nbinomial", offsetScale=1000^2, timeModel="rw2"))
  model$estimate()
  model$collectEstimates()
  model$collectHyperparameters()
  summary(model$result)
  model$getEstimatesFileName()
  
  #model <- FinlandSmoothModelTemporal$new(study=study)
  #model$setModelName("nbinomial", timeModels[task_id])
  #model$loadEstimates()
  
  habitatWeightsRaster <- study$loadHabitatWeightsRaster()
  model$collectEstimates()
  model$collectHyperparameters()
  populationDensity <- model$getPopulationDensity(templateRaster=habitatWeightsRaster, getSD=FALSE)
  populationDensity$mean$weight(habitatWeightsRaster)
  populationSize <- populationDensity$mean$integrate(volume=FinlandPopulationSize$new(study=study))
  populationSize$loadValidationData()
  populationSize
} else if (FALSE) {
  # OLD
  
  # For the full estimation
  
  estimate <- function(response, spatialModel, timeModel) {
    context <- Context$new(resultDataDirectory=wd.data.results, processedDataDirectory=wd.data.processed, rawDataDirectory=wd.data.raw, scratchDirectory=wd.scratch, figuresDirectory=wd.figures)
    study <- FinlandWTCStudy$new(context=context, response=response, distanceCovariatesModel=~populationDensity+rrday+snow+tday-1, trackSampleInterval=2)
    
    model <- if (spatialModel) FinlandSmoothModelSpatioTemporal$new(study=study)
    else FinlandSmoothModelTemporal$new(study=study)
    
    meshParams <- list(coordsScale=1e-6, maxEdge=c(.1e6, .2e6), cutOff=.01e6)
    params <- if (timeModel == "ar2") list(family="nbinomial", offsetScale=1000^2, meshParams=meshParams, timeModel="ar", model=response ~ 1 + f(year, model="ar", order=2))
    else list(family="nbinomial", offsetScale=1000^2, meshParams=meshParams, timeModel=timeModel)
    
    model <- study$estimate(model=model, params=params, save=T)
  }
  
  responses <- c("canis.lupus", "lynx.lynx", "rangifer.tarandus.fennicus")
  response <- responses[task_id]
  spatialModels <- c(T, T, F)
  timeModels <- c("ar1", "ar1", "rw2")
  estimate(response=response, spatialModel=spatialModels[task_id], timeModel=timeModels[task_id])
} else {
  estimate <- function(response, modelName) {
    context <- Context(resultDataDirectory=wd.data.results, processedDataDirectory=wd.data.processed, rawDataDirectory=wd.data.raw, scratchDirectory=wd.scratch, figuresDirectory=wd.figures)
    study <- FinlandWTCStudy(context=context, response=response, distanceCovariatesModel=~populationDensity+rrday+snow+tday-1, trackSampleInterval=2)
    
    if (modelName == "SmoothModel-nbinomial-matern-ar1") {
      model <- FinlandSmoothModelSpatioTemporal(study=study)
      modelParams <- list(family="nbinomial", offsetScale=1000^2, meshParams=study$studyArea$getMesh(), timeModel="ar1")      
    }
    else if (modelName == "SmoothModel-nbinomial-ar1") {
      model <- FinlandSmoothModelTemporal(study=study)
      modelParams <- list(family="nbinomial", offsetScale=1000^2, timeModel="ar1")
    }
    else if (modelName == "FMPModel") {
      model <- FinlandFMPModel(study=study)
      modelParams <- NULL
    }
    
    model <- study$estimate(model=model, params=modelParams, save=T)
  }
  
  responses <- c("canis.lupus", "lynx.lynx", "rangifer.tarandus.fennicus")
  response <- responses[task_id]
  #spatialModels <- c(T, T, F)
  #timeModels <- c("ar1", "ar1", "rw2")
  estimate(response=response, modelName)
}
