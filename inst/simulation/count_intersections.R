# Run test:
# ./parallel_r.py -t 1:5 -n 3 -l 10.0 -b ~/tmp/blacklist.txt -v ~/git/Winter-Track-Counts/inst/simulation/count_intersections.R test A
# Run full:
# ./parallel_r.py -t 1:50 -n 50 -l 10.0 -b ~/tmp/blacklist.txt -v ~/git/Winter-Track-Counts/inst/simulation/count_intersections.R notest A

# library(devtools); install_github("statguy/Winter-Track-Counts")


countIntersections <- function(scenario, iteration, test) {
  context <- Context$new(resultDataDirectory=wd.data.results, processedDataDirectory=wd.data.processed, rawDataDirectory=wd.data.raw, scratchDirectory=wd.scratch, figuresDirectory=wd.figures)
  mss <- if (scenario == "A") MovementSimulationScenarioA$new()$newInstance(context=context) else stop("Unknown scenario ", scenario)
  study <- mss$study
  tracks <- study$loadTracks(iteration=iteration)
  surveyRoutes <- study$loadSurveyRoutes()
  intersections <- SimulatedIntersections$new(study=study, iteration=iteration)
  if (!test) {
    intersections$findIntersections(tracks, surveyRoutes, dimension=1)
    intersections$saveIntersections()
  }
  else {
    surveyRoutes$surveyRoutes <- surveyRoutes$surveyRoutes[1:2]
    surveyRoutes$centroids <- surveyRoutes$centroids[1:2]
    surveyRoutes$lengths <- surveyRoutes$lengths[1:2]
    intersections$findIntersections(tracks, surveyRoutes, dimension=1)
    message("SUCCESS")
  }
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 3) stop("Invalid arguments.")
test <- args[1]
scenario <- args[2]
task_id <- args[length(args)]
message("Arguments provided:")
print(args)

library(parallel)
library(doMC)
registerDoMC(cores=detectCores())

library(WTC)
source("~/git/Winter-Track-Counts/setup/WTC-Boot.R")

countIntersections(scenario=scenario, iteration=as.integer(task_id), test=test=="test")

traceback()
