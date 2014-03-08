# Run test:
# ./parallel_r.py -t 1:5 -n 2 -l 10.0 -b ~/tmp/blacklist.txt -v ~/git/Winter-Track-Counts/inst/simulation/count_intersections.R test
# Run full:
# ./parallel_r.py -t 1:50 -n 50 -l 10.0 -b ~/tmp/blacklist.txt -v ~/git/Winter-Track-Counts/inst/simulation/count_intersections.R notest

err <- try({
  countIntersections <- function(iteration, test) {
    context <- Context$new(resultDataDirectory=wd.data.results, processedDataDirectory=wd.data.processed, rawDataDirectory=wd.data.raw, scratchDirectory=wd.scratch, figuresDirectory=wd.figures)
    mss <- MovementSimulationScenarioA$new()$newInstance(context=context)
    study <- mss$study
    tracks <- study$loadTracks(iteration=iteration)
    surveyRoutes <- study$loadSurveyRoutes()
    intersections <- SimulatedIntersections$new(study=study, iteration=iteration)
    if (!test) {
      intersections$findIntersections(tracks, surveyRoutes, dimension=1)
      intersections$saveIntersections()
    }
    else message("SUCCESS")
  }
  
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
  
  countIntersections(iteration=task_id, test=test=="test")
})

if (inherits(err, "try-error")) {
  traceback()
  message(err); stop("count_intersections.R; err = 1, msg = ", err[1])
}
