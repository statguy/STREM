# Run test:
# ./parallel_r.py -t 1:2 -n 2 -l 10.0 -b ~/tmp/blacklist.txt -v ~/git/Winter-Track-Counts/inst/simulation/simulate.R test A
# Run full:
# ./parallel_r.py -t 1:50 -n 50 -l 10.0 -b ~/tmp/blacklist.txt -v ~/git/Winter-Track-Counts/inst/simulation/simulate.R notest A

# library(devtools); install_github("statguy/Winter-Track-Counts")

simulate = function(scenario, iteration, nAgents=as.integer(200), nYears=as.integer(20),
                    nDays=as.integer(365),
                    CRWCorrelation=0.7,
                    BCRWCorrelationBiasTradeoff=0.7,
                    isTest=FALSE,
                    plot=FALSE,
                    returnMSS=FALSE) {
  library(WTC)
  source("~/git/Winter-Track-Counts/setup/WTC-Boot.R")
  
  context <- Context$new(resultDataDirectory=wd.data.results, processedDataDirectory=wd.data.processed, rawDataDirectory=wd.data.raw, scratchDirectory=wd.scratch, figuresDirectory=wd.figures)
  mss <- {
    if (scenario == "A") MovementSimulationScenarioA$new(nAgents=nAgents, years=nYears, days=nDays, CRWCorrelation=CRWCorrelation)$setup(context=context, isTest=isTest)
    else if (scenario == "B") MovementSimulationScenarioB$new(nAgents=nAgents, years=nYears, days=nDays, CRWCorrelation=CRWCorrelation, BCRWCorrelationBiasTradeoff=BCRWCorrelationBiasTradeoff)$setup(context=context, isTest=isTest)
    else if (scenario == "C") MovementSimulationScenarioC$new(nAgents=as.integer(nAgents/5), years=nYears, days=nDays, CRWCorrelation=CRWCorrelation)$setup(context=context, isTest=isTest)
    else if (scenario == "D") MovementSimulationScenarioD$new(nAgents=nAgents, years=nYears, days=nDays, CRWCorrelation=CRWCorrelation)$setup(context=context, isTest=isTest)
    else if (scenario == "E") MovementSimulationScenarioE$new(nAgents=nAgents, years=nYears, days=nDays, CRWCorrelation=CRWCorrelation)$setup(context=context, readHabitatIntoMemory=T, isTest=isTest)
    else if (scenario == "F") MovementSimulationScenarioF$new(nAgents=as.integer(nAgents/5), years=nYears, days=nDays, CRWCorrelation=CRWCorrelation, BCRWCorrelationBiasTradeoff=BCRWCorrelationBiasTradeoff)$setup(context=context, readHabitatIntoMemory=T, isTest=isTest)
    else stop("unsupported")
  }

  if (returnMSS) return(mss)
    
  study <- mss$study
  #message("Study area = ", study$studyArea$boundary@polygons[[1]]@area / 1000^2, " km^2")
  tracks <- mss$simulate(iteration=iteration, save=T)
}


args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 3) stop("Invalid arguments.")
isTest <- args[1] == "test"
scenario <- args[2]
task_id <- args[length(args)]
message("Arguments provided:")
print(args)

library(parallel)
library(doMC)
registerDoMC(cores=detectCores())
library(WTC)
source("~/git/Winter-Track-Counts/setup/WTC-Boot.R")

#for (task_id in 1:100)
if (isTest) simulate(scenario=scenario, iteration=as.integer(task_id), nAgents=as.integer(50), nYears=as.integer(1), nDays=as.integer(59), isTest=TRUE)
else simulate(scenario=scenario, iteration=as.integer(task_id))
