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

# For testing

if (test == "test") {
  context <- Context$new(resultDataDirectory=wd.data.results, processedDataDirectory=wd.data.processed, rawDataDirectory=wd.data.raw, scratchDirectory=wd.scratch, figuresDirectory=wd.figures)
  study <- FinlandWTCStudy$new(context=context, response="lynx.lynx")
  model <- study$estimate(quick=TRUE)
  surveyRoutes <- study$loadSurveyRoutes()
  model$plotMesh(surveyRoutes)
  model <- study$estimate(test=TRUE)
  model <- study$estimate(quick=TRUE)
}
else {
# For full estimation
  estimate <- function(response) {
    context <- Context$new(resultDataDirectory=wd.data.results, processedDataDirectory=wd.data.processed, rawDataDirectory=wd.data.raw, scratchDirectory=wd.scratch, figuresDirectory=wd.figures)
    study <- FinlandWTCStudy$new(context=context, response=response)
    model <- study$estimate()      
  }
  
  response <- if (task_id == 1) "canis.lupus"
  else if (task_id == 2) "lynx.lynx"
  else if (task_id == 3) "rangifer.tarandus.fennicus"
}
