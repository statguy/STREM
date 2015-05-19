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
  context <- Context(resultDataDirectory=wd.data.results, processedDataDirectory=wd.data.processed, rawDataDirectory=wd.data.raw, scratchDirectory=wd.scratch, figuresDirectory=wd.figures)
  study <- FinlandWTCStudy(context=context, response=response)
  
  intersections <- FinlandWTCIntersections$new(study=study)$setCovariatesId(tag="density")$loadCovariates()$removeMissingCovariates()
  covariates <- intersections$compressCovariates()
  #cor(covariates[,fitCovariates]) > .5
  fitCovariates <- stringr::str_subset(colnames(covariates), "^PC\\d+_habitat_+|^topography_rel_+")
  #covariatesModel <- reformulate(fitCovariates, intercept=FALSE)
  covariatesModel <- ~ PC1_habitat_62.5 + PC2_habitat_62.5 + topography_rel_1000
  intersections$intersections@data <- cbind(intersections$getData(), covariates[,fitCovariates])
  model <- SmoothModelSpatioTemporal$new(study=study)
  modelParams <- study$getModelParams(modelName="SmoothModel-nbinomial-matern-ar1")
  modelParams$meshParams$cutOff <- 600000
  modelParams$maxEdge <- c(4e5,5e5)
  model$setup(intersections=intersections, params=modelParams, covariatesModel=covariatesModel, tag="full")
  model$estimate()
  summary(model$result)
  
  
  covariatesModel <- ~ PC1_habitat_1000 + PC2_habitat_1000 + topography_rel_4000
  model$setup(intersections=intersections, params=modelParams, covariatesModel=covariatesModel, tag="full")
  model$estimate()
  summary(model$result)
 
  
} else {
  estimate <- function(response, modelName) {
    context <- Context(resultDataDirectory=wd.data.results, processedDataDirectory=wd.data.processed, rawDataDirectory=wd.data.raw, scratchDirectory=wd.scratch, figuresDirectory=wd.figures)
    study <- FinlandWTCStudy(context=context, response=response)
    model <- study$getModel(modelName=modelName)
    modelParams <- study$getModelParams(modelName=modelName)    
    model <- study$estimate(model=model, params=modelParams, save=T)
  }
  
  responses <- c("canis.lupus", "lynx.lynx", "rangifer.tarandus.fennicus")
  response <- responses[task_id]
  estimate(response=response, modelName=modelName)
}
