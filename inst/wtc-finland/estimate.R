# ./parallel_r.py -t 1:3 -n 6 -l 10.0 -b ~/tmp/blacklist.txt -v ~/git/Winter-Track-Counts/inst/wtc/estimate.R notest FMPModel
# ./parallel_r.py -t 3 -n 2 -l 10.0 -b ~/tmp/blacklist.txt -v ~/git/Winter-Track-Counts/inst/wtc/estimate.R notest FMPModel

# library(devtools); install_github("statguy/Winter-Track-Counts")

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 3) stop("Invalid arguments.")
test <- args[1]
modelName <- args[2]
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

  context <- Context(resultDataDirectory=wd.data.results, processedDataDirectory=wd.data.processed, rawDataDirectory=wd.data.raw, scratchDirectory=wd.scratch, figuresDirectory=wd.figures)
  study <- FinlandWTCStudy(context=context, response=response, distanceCovariatesModel=~populationDensity+rrday+snow+tday-1, trackSampleInterval=2)
  
  intersections <- study$loadIntersections()
  model <- FinlandFMPModel(study=study)
  model$setup(intersections=intersections, params=NULL)
  model$estimate()
  model$collectEstimates()
  #model$collectHyperparameters()
  #summary(model$result)
  model$getEstimatesFileName()
  
  #model <- FinlandSmoothModelTemporal$new(study=study)
  #model$setModelName("nbinomial", timeModels[task_id])
  #model$loadEstimates()
  
  model$collectEstimates()
  #model$collectHyperparameters()
  populationDensity <- model$getPopulationDensity(templateRaster=habitatWeightsRaster, getSD=FALSE)
  populationSize <- FinlandPopulationSize(study=study, modelName=modelName)$getPopulationSize(populationDensity=populationDensity$mean)
  populationSize
} else {
  estimate <- function(response, modelName) {
    context <- Context(resultDataDirectory=wd.data.results, processedDataDirectory=wd.data.processed, rawDataDirectory=wd.data.raw, scratchDirectory=wd.scratch, figuresDirectory=wd.figures)
    study <- FinlandWTCStudy(context=context, response=response, distanceCovariatesModel=~populationDensity+rrday+snow+tday-1, trackSampleInterval=2)
    model <- study$getModel(modelName=modelName)
    modelParams <- study$getModelParams(modelName=modelName)    
    model <- study$estimate(model=model, params=modelParams, save=T)
  }
  
  responses <- c("canis.lupus", "lynx.lynx", "rangifer.tarandus.fennicus")
  response <- responses[task_id]
  estimate(response=response, modelName)
}
