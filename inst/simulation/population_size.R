# Run full:
# ./parallel_r.py -t 1:50 -n 60 -l 10.0 -b ~/tmp/blacklist.txt -v ~/git/Winter-Track-Counts/inst/simulation/population_size.R notest A
# R --vanilla --args notest A 1 < population_size.R

# library(devtools); install_github("statguy/Winter-Track-Counts")
# echo 'library(devtools); install_github("statguy/Winter-Track-Counts")' | R --slave
# args <- c("notest","A","1")


population_size <- function(mss, iteration, test) {
  if (test) {
    study <- mss$study
    #iteration <- as.integer(50)
    estimates <- study$loadEstimates(iteration=iteration)
    estimates$collectEstimates(predictionWeights=estimates$getPredictedOffset())
    #estimates$collectHyperparameters()
    populationSize <- estimates$getPopulationSize(withHabitatWeights=mss$hasHabitatWeights())
    populationSize$loadValidationData()
    populationSize
    colSums(populationSize$sizeData[,-1])
    colMeans(populationSize$sizeData[,-1])
    populationSize$plotPopulationSize()
  }
  else {
    study <- mss$study
    estimates <- SimulatedSmoothModelSpatioTemporal(study=study, iteration=iteration)
    estimates$setModelName(family="nbinomial", randomEffect=paste("matern", "ar1", sep="-"))
    study$getPopulationSize(estimates=estimates)
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


context <- Context$new(resultDataDirectory=wd.data.results, processedDataDirectory=wd.data.processed, rawDataDirectory=wd.data.raw, scratchDirectory=wd.scratch, figuresDirectory=wd.figures)
mss <- {
  if (scenario == "A") MovementSimulationScenarioA$new()$setup(context=context)
  else if (scenario == "B") MovementSimulationScenarioB$new()$setup(context=context)
  else if (scenario == "C") MovementSimulationScenarioC$new()$setup(context=context)
  else if (scenario == "D") MovementSimulationScenarioD$new()$setup(context=context)
  else if (scenario == "E") MovementSimulationScenarioE$new()$setup(context=context)
  else if (scenario == "F") MovementSimulationScenarioF$new()$setup(context=context)
  else stop("unsupported")
}

population_size(mss=mss, iteration=as.integer(task_id), test=test=="test")
