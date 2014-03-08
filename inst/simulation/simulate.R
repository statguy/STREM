# Run test:
# ./parallel_r.py -t 1:5 -n 2 -l 10.0 -b ~/tmp/blacklist.txt -v ~/git/Winter-Track-Counts/inst/simulation/simulate.R test test
# Run full:
# ./parallel_r.py -t 1:50 -n 50 -l 10.0 -b ~/tmp/blacklist.txt -v ~/git/Winter-Track-Counts/inst/simulation/simulate.R notest A

dryRun <- function(test) {
  context <- Context$new(resultDataDirectory=wd.data.results, processedDataDirectory=wd.data.processed, rawDataDirectory=wd.data.raw, scratchDirectory=wd.scratch, figuresDirectory=wd.figures)
  mss <- if (test == "test") MovementSimulationScenarioA$new(nAgents=as.integer(2), nIterations=as.integer(5), years=as.integer(2))$newInstance(context=context)
}

simulateA <- function(test)
  context <- Context$new(resultDataDirectory=wd.data.results, processedDataDirectory=wd.data.processed, rawDataDirectory=wd.data.raw, scratchDirectory=wd.scratch, figuresDirectory=wd.figures)
  mss <- if (test == "test") MovementSimulationScenarioA$new(nAgents=as.integer(2), nIterations=as.integer(5), years=as.integer(2))$newInstance(context=context)
  else MovementSimulationScenarioA$new()$newInstance(context=context)
  mss$simulateSingle(iteration=task_id)
}

err <- try({
  args <- commandArgs(trailingOnly=TRUE)
  if (length(args) != 3) stop("Invalid arguments.")
  test <- args[1]
  scenario <- args[2]
  task_id <- args[length(args)]

  library(parallel)
  library(doMC)
  registerDoMC(cores=detectCores())

  library(WTC)
  source("~/git/Winter-Track-Counts/setup/WTC-Boot.R")

  if (scenario == "test") dryRun(test=test)
  else if (scenario == "A") simulateA(test=test)
  else stop("Unsupported scenario ", scenario)
})

if (inherits(err, "try-error")) {
  traceback()
  message(err); stop("simulate.R; err = 1, msg = ", err[1])
}
