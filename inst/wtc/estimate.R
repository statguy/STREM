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
  study <- FinlandWTCStudy$new(context=context, response="lynx.lynx", distanceCovariatesModel=~populationDensity+rrday+snow+tday-1, trackSampleInterval=2)
  
  intersections <- study$loadIntersections()
  model <- FinlandSmoothModel$new(study = study)
  model$setup(intersections = intersections, meshParams = meshParams, useCovariates=F)
  
  
  meshParams <- list(maxEdge=c(.2e6, .4e6), cutOff=.1e6, coordsScale=1e-6)
  #meshParams <- list(maxEdge=c(.06e6, .15e6), cutOff=.03e6, coordsScale=1e-6)
  model <- study$estimate(meshParams=meshParams, interceptPriorParams=list(mean=200, sd=199), test=TRUE)
  surveyRoutes <- study$loadSurveyRoutes()
  model$plotMesh(surveyRoutes)
  model <- study$estimate(meshParams=meshParams, interceptPriorParams=list(mean=200, sd=199))
  
  intersections <- study$loadIntersections()
  intersections$predictDistances()
  
  model$collectEstimates()
  
} else {
  # For full estimation
  
  estimate <- function(response) {
    context <- Context$new(resultDataDirectory=wd.data.results, processedDataDirectory=wd.data.processed, rawDataDirectory=wd.data.raw, scratchDirectory=wd.scratch, figuresDirectory=wd.figures)
    study <- FinlandWTCStudy$new(context=context, response=response, distanceCovariatesModel=~populationDensity+rrday+snow+tday-1, trackSampleInterval=2)
    
    #meshParams <- list(maxEdge=c(.07e6, .2e6), cutOff=.01e6, coordsScale=1e-6)
    meshParams <- list(maxEdge=c(.1e6, .2e6), cutOff=.01e6, coordsScale=1e-6)
    model <- study$estimate(meshParams=meshParams, interceptPriorParams=list(mean=200, sd=199), save=TRUE)
  }
  
  response <- if (task_id == 1) "canis.lupus"
  else if (task_id == 2) "lynx.lynx"
  else if (task_id == 3) "rangifer.tarandus.fennicus"
  estimate(response=response)
}
