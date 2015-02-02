# Run test:
# ./parallel_r.py -t 1:50 -n 60 -l 10.0 -b ~/tmp/blacklist.txt -v ~/git/Winter-Track-Counts/inst/simulation/count_intersections.R test A
# Run full:
# ./parallel_r.py -t 1:50 -n 60 -l 10.0 -b ~/tmp/blacklist.txt -v ~/git/Winter-Track-Counts/inst/simulation/count_intersections.R notest A
# ./parallel_r.py -t 1:5 -n 6 -l 10.0 -b ~/tmp/blacklist.txt -v ~/git/Winter-Track-Counts/inst/simulation/count_intersections.R notest Acombined

# library(devtools); install_github("statguy/Winter-Track-Counts")

countIntersections <- function(scenario, iteration, nSurveyRoutes=as.integer(500), countDays=as.integer(1), isTest=F) {
  mss <- getMSS(scenario=scenario, nSurveyRoutes=nSurveyRoutes, isTest=isTest)
  study <- mss$study
  study$countIntersections(surveyRoutes=mss$getSurveyRoutes(), iteration=iteration, days=countDays)
}

library(parallel)
library(doMC)
registerDoMC(cores=detectCores())
library(WTC)
source("~/git/Winter-Track-Counts/setup/WTC-Boot.R")

parseArguments()
countDays <- if (!is.null(extraArgs)) as.integer(extraArgs[1]) else as.integer(1)

message("Counting intersections for scenario = ", scenario, " iteration = ", task_id, " days = ", countDays)
{
  if (isTest) countIntersections(scenario=scenario, iteration=as.integer(task_id), nSurveyRoutes=as.integer(50), countDays=countDays, isTest=T)
  else countIntersections(scenario=scenario, iteration=as.integer(task_id), countDays=countDays, isTest=F)
}
