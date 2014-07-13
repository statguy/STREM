# Run test:
# ./parallel_r.py -t 1:2 -n 2 -l 10.0 -b ~/tmp/blacklist.txt -v ~/git/Winter-Track-Counts/inst/simulation/simulate.R test A
# Run full:
# ./parallel_r.py -t 1:50 -n 50 -l 10.0 -b ~/tmp/blacklist.txt -v ~/git/Winter-Track-Counts/inst/simulation/simulate.R notest A


simulate = function(scenario, iteration, nSurveyRoutes=as.integer(50), nAgents=as.integer(20), nYears=as.integer(20),
                    nDays=as.integer(365), nIterations=as.integer(100),
                    meshParams=list(coordsScale=1e-6, maxEdge=c(.01e6, .02e6), cutOff=.007e6),
                    modelParams=list(family="nbinomial", offsetScale=1000^2, meshParams=meshParams, timeModel="ar1"),
                    CRWCorrelation=0.7,
                    BCRWCorrelationBiasTradeoff=0.7,
                    plot=FALSE) {
  library(WTC)
  source("~/git/Winter-Track-Counts/setup/WTC-Boot.R")
  
  context <- Context$new(resultDataDirectory=wd.data.results, processedDataDirectory=wd.data.processed, rawDataDirectory=wd.data.raw, scratchDirectory=wd.scratch, figuresDirectory=wd.figures)
  mss <- {
    if (scenario == "A") MovementSimulationScenarioA$new(nAgents=nAgents, years=nYears, days=nDays, CRWCorrelation=CRWCorrelation)$setup(context=context)
    else if (scenario == "B") MovementSimulationScenarioB$new(nAgents=nAgents, years=nYears, days=nDays, CRWCorrelation=CRWCorrelation, BCRWCorrelationBiasTradeoff=BCRWCorrelationBiasTradeoff)$setup(context=context)
    else if (scenario == "C") MovementSimulationScenarioC$new(nAgents=nAgents, years=nYears, days=nDays, CRWCorrelation=CRWCorrelation)$setup(context=context)
    else if (scenario == "D") MovementSimulationScenarioD$new(nAgents=nAgents, years=nYears, days=nDays, CRWCorrelation=CRWCorrelation)$setup(context=context)
    else if (scenario == "E") MovementSimulationScenarioE$new(nAgents=nAgents, years=nYears, days=nDays, CRWCorrelation=CRWCorrelation, nSurveyRoutes=nSurveyRoutes)$setup(context=context)
    else if (scenario == "F") MovementSimulationScenarioF$new(nAgents=nAgents, years=nYears, days=nDays, CRWCorrelation=CRWCorrelation, BCRWCorrelationBiasTradeoff=BCRWCorrelationBiasTradeoff, nSurveyRoutes=nSurveyRoutes)$setup(context=context)
    else stop("unsupported")
  }
  
  study <- mss$study
  #message("Study area = ", study$studyArea$boundary@polygons[[1]]@area / 1000^2, " km^2")
  tracks <- mss$simulate(iteration=iteration, save=T)
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

simulate(scenario=scenario, iteration=as.integer(task_id))
