# Run test:
# ./parallel_r.py -t 1:3 -n 3 -l 10.0 -b ~/tmp/blacklist.txt -v ~/git/Winter-Track-Counts/inst/simulation/count_intersections.R test A
# Run full:
# ./parallel_r.py -t 1:50 -n 50 -l 10.0 -b ~/tmp/blacklist.txt -v ~/git/Winter-Track-Counts/inst/simulation/count_intersections.R notest A

# library(devtools); install_github("statguy/Winter-Track-Counts")


countIntersections <- function(study, iteration, test) {
  tracks <- study$loadTracks(iteration=iteration)
  if (test) {
    surveyRoutes <- study$loadSurveyRoutes()
    intersections <- SimulatedIntersections$new(study=study, iteration=iteration)
    surveyRoutes$surveyRoutes <- surveyRoutes$surveyRoutes[1:3]
    surveyRoutes$centroids <- surveyRoutes$centroids[1:3]
    surveyRoutes$lengths <- surveyRoutes$lengths[1:3]
    observationTracks <- tracks$randomizeObservationDayTracks()
    intersections$findIntersections(observationTracks, surveyRoutes, dimension=1)
    message("SUCCESS")
  }
  else tracks$countIntersections()
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

context <- Context$new(resultDataDirectory=wd.data.results, processedDataDirectory=wd.data.processed, rawDataDirectory=wd.data.raw, scratchDirectory=wd.scratch, figuresDirectory=wd.figures)
mss <- if (scenario == "A") MovementSimulationScenarioA$new()$newInstance(context=context) else stop("Unknown scenario ", scenario)
study <- mss$study
countIntersections(study=study, iteration=as.integer(task_id), test=test=="test")
