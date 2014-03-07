# Run test:
# ./parallel_ssh_r.py -n 5 -m 2 -l 10.0 -b ~/tmp/blacklist.txt -v ~/git/Winter-Track-Counts/inst/simulation/simulation.R test A
# Run full:
# ./parallel_ssh_r.py -n 50 -m 50 -l 10.0 -b ~/tmp/blacklist.txt -v ~/git/Winter-Track-Counts/inst/simulation/simulation.R notest A

function <- simulateA(test)
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

  if (scenario == "A") simulateA(test=test)
  else stop("Unsupported scenario ", scenario)
})

if (inherits(err, "try-error")) {
   traceback()
   message(err); stop("simulation.R; err = 1, msg = ", err[1])
}
