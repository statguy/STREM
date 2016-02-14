# ./parallel_r.py -t 1:3 -n 6 -l 10.0 -b ~/tmp/blacklist.txt -v ~/git/Winter-Track-Counts/inst/wtc/estimate.R notest FMPModel
# ./parallel_r.py -t 3 -n 2 -l 10.0 -b ~/tmp/blacklist.txt -v ~/git/Winter-Track-Counts/inst/wtc/estimate.R notest FMPModel

# library(devtools); install_github("statguy/WTC")

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
library(STREM)
source("~/git/STREM/setup/WTC-Boot.R")

if (test == "test") {
  response <- "canis.lupus"
  
  context <- Context(resultDataDirectory=wd.data.results, processedDataDirectory=wd.data.processed, rawDataDirectory=wd.data.raw, scratchDirectory=wd.scratch, figuresDirectory=wd.figures)
  study <- FinlandWTCStudy(context=context, response=response)
  
  intersections <- FinlandWTCIntersections$new(study=study)$setCovariatesId(tag="density")$loadCovariates()$removeMissingCovariates()
  #plot(study$studyArea$boundary)
  #lines(matrix(c(3186635,7300000,3732935,7300000), byrow=T, ncol=2))
  index <- coordinates(intersections$intersections)[,2] < 7.3e6 # Cut off Lapland since all wolves are killed there
  index <- index & intersections$intersections$duration == 1 # Use durations of 1 only
  intersections$intersections <- intersections$intersections[index,]
  colnames(intersections$intersections@data) <- stringr::str_replace_all(colnames(intersections$intersections@data), "[\\[\\]_]", ".") # there is a problem with column names in INLA or somewhere
  intersections$covariateNames <- stringr::str_replace_all(intersections$covariateNames, "[\\[\\]_]", ".")

  #ddply(intersections$intersections@data, .(year), nrow)
  
  dataNames <- colnames(intersections$intersections@data)
  for (scale in 2^(0:10) * 62.5) {
    #variables <- stringr::str_subset(dataNames, "^habitat\\.\\w+\\.\\.62\\.5$")
    #intersections$compressCovariates(variables, name="habitat.62.5")
    variables <- stringr::str_subset(dataNames, paste0("^habitat\\.\\w+\\.\\.", scale, "$"))
    intersections$compressCovariates(variables, name=paste0("habitat.", scale))
  }
  variables <- stringr::str_subset(dataNames, "^topography\\.rel\\.\\d+$")
  intersections$scaleCovariates(variables)
  
#  model <- SmoothModelSpatioTemporal$new(study=study)
#  modelParams <- study$getModelParams(modelName="SmoothModel-nbinomial-matern-ar1")
#  #modelParams$meshParams$cutOff <- 300000
#  #modelParams$maxEdge <- c(3e5,5e5)
#  covariatesModel <- ~ habitat.Agriculture..8000 + habitat.Urban..250 + habitat.Forestland..32000 + topography.rel.32000
#  #covariatesModel <- ~ topography_rel_64000
#  model$setup(intersections=intersections, params=modelParams, covariatesModel=covariatesModel, tag="full", timeout=120)
#  print(model$model)
#  model$estimate()
#  summary(model$result)
  
#  model$setup(intersections=intersections, params=modelParams, tag="full")
#  print(model$model)
#  model$estimate()
#  summary(model$result)
 
  
  reload(inst("SpaceTimeModels"))
  library(SpaceTimeModels)

  data <- ddply(cbind(coordinates(intersections$intersections), intersections$intersections@data), .(id),
                function(x) data.frame(x=x$coords.x1[1], y=x$coords.x2[1], length=sum(x$length), distance=sum(x$distance), intersections=sum(x$intersections), x[1,10:ncol(x)]))
  
  mesh <- NonConvexHullMesh$new(SpatialPoints(coordinates(intersections$intersections)), knotsScale=1e6, cutoff=0.01e6, offset=c(1,0.1)*1e6, maxEdge=c(0.05,1)*1e6, convex=0.12)
  mesh$getINLAMesh()$n
  
  meanDistance <- mean(intersections$intersections@data$distance)
  
  model <- ContinuousSpaceModel$new()
  model$setSpatialMesh(mesh)
  model$setSpatialPrior()
  model$setSmoothingModel()
  #model$addObservationStack(sp=SpatialPoints(data[,c("x","y")]), response=data$intersections, offset=data$length * data$distance, covariates=data, tag="obs")
  model$addObservationStack(sp=SpatialPoints(data[,c("x","y")]), response=data$intersections, offset=data$length * meanDistance, covariates=data, tag="obs")
  model$setLikelihood("nbinomial")
  model$estimate(verbose=F)
  summary(model$getResult())

  names(data)
  model$setCovariatesModel(~ habitat.Agriculture..8000 + habitat.Urban..250 + habitat.Forestland..32000 + topography.rel.32000, data)
  #model$clearStack()$addObservationStack(sp=SpatialPoints(data[,c("x","y")]), response=data$intersections, offset=data$length * data$distance, covariates=data, tag="obs")
  model$clearStack()$addObservationStack(sp=SpatialPoints(data[,c("x","y")]), response=data$intersections, offset=data$length * meanDistance, covariates=data, tag="obs")
  model$estimate(verbose=F)
  summary(model$getResult())
  
  
  
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
