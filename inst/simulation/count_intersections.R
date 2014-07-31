# Run test:
# ./parallel_r.py -t 1:3 -n 3 -l 10.0 -b ~/tmp/blacklist.txt -v ~/git/Winter-Track-Counts/inst/simulation/count_intersections.R test A
# Run full:
# ./parallel_r.py -t 1:50 -n 60 -l 10.0 -b ~/tmp/blacklist.txt -v ~/git/Winter-Track-Counts/inst/simulation/count_intersections.R notest A

# library(devtools); install_github("statguy/Winter-Track-Counts")

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 3) stop("Invalid arguments.")
test <- args[1] == "test"
scenario <- args[2]
task_id <- args[length(args)]
message("Arguments provided:")
print(args)

library(parallel)
library(doMC)
registerDoMC(cores=detectCores())

library(WTC)
source("~/git/Winter-Track-Counts/setup/WTC-Boot.R")


countIntersections <- function(scenario, iteration, nSurveyRoutes=as.integer(500), isTest=F) {
  context <- Context$new(resultDataDirectory=wd.data.results, processedDataDirectory=wd.data.processed, rawDataDirectory=wd.data.raw, scratchDirectory=wd.scratch, figuresDirectory=wd.figures)
  mss <- {
    if (scenario == "A") MovementSimulationScenarioA$new()$setup(context=context, nSurveyRoutes=nSurveyRoutes, isTest=isTest)
    else if (scenario == "B") MovementSimulationScenarioB$new()$setup(context=context, nSurveyRoutes=nSurveyRoutes, isTest=isTest)
    else if (scenario == "C") MovementSimulationScenarioC$new()$setup(context=context, nSurveyRoutes=nSurveyRoutes, isTest=isTest)
    else if (scenario == "D") MovementSimulationScenarioD$new()$setup(context=context, nSurveyRoutes=nSurveyRoutes, isTest=isTest)
    else if (scenario == "E") MovementSimulationScenarioE$new()$setup(context=context, nSurveyRoutes=nSurveyRoutes, isTest=isTest)
    else if (scenario == "F") MovementSimulationScenarioF$new()$setup(context=context, nSurveyRoutes=nSurveyRoutes, isTest=isTest)
    else stop("unsupported")
  }
  
  study <- mss$study
  study$countIntersections(surveyRoutes=mss$getSurveyRoutes(), iteration=iteration)
}

#scenario<-"A"; task_id<-1; test<-F

if (test) countIntersections(scenario=scenario, iteration=as.integer(task_id), nSurveyRoutes=as.integer(50), isTest=T)
else countIntersections(scenario=scenario, iteration=as.integer(task_id), isTest=F)
