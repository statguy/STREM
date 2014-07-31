# Run test:
# ./parallel_r.py -t 1:50 -n 60 -l 10.0 -b ~/tmp/blacklist.txt -v ~/git/Winter-Track-Counts/inst/simulation/population_size.R test A
# Run full:
# ./parallel_r.py -t 1:50 -n 60 -l 10.0 -b ~/tmp/blacklist.txt -v ~/git/Winter-Track-Counts/inst/simulation/population_size.R notest A
# R --vanilla --args notest A 1 < population_size.R

# library(devtools); install_github("statguy/Winter-Track-Counts")
# echo 'library(devtools); install_github("statguy/Winter-Track-Counts")' | R --slave
# args <- c("notest","A","1")


population_size <- function(mss, iteration, isTest, otherTest=F) {
  context <- Context$new(resultDataDirectory=wd.data.results, processedDataDirectory=wd.data.processed, rawDataDirectory=wd.data.raw, scratchDirectory=wd.scratch, figuresDirectory=wd.figures)
  mss <- {
    if (scenario == "A") MovementSimulationScenarioA$new()$setup(context=context, isTest=isTest)
    else if (scenario == "B") MovementSimulationScenarioB$new()$setup(context=context, isTest=isTest)
    else if (scenario == "C") MovementSimulationScenarioC$new()$setup(context=context, isTest=isTest)
    else if (scenario == "D") MovementSimulationScenarioD$new()$setup(context=context, isTest=isTest)
    else if (scenario == "E") MovementSimulationScenarioE$new()$setup(context=context, isTest=isTest)
    else if (scenario == "F") MovementSimulationScenarioF$new()$setup(context=context, isTest=isTest)
    else stop("unsupported")
  }
  study <- mss$study
  
  if (otherTest) {
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
    estimates <- SimulatedSmoothModelSpatioTemporal(study=study, iteration=iteration)
    estimates$setModelName(family="nbinomial", randomEffect=paste("matern", "ar1", sep="-"))
    habitatWeights <- study$getHabitatWeights(iteration=iteration)
    populationSize <- study$getPopulationSize(estimates=estimates, habitatWeights=habitatWeights)
    return(invisible(populationSize))
  }
}

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

x <- population_size(mss=mss, iteration=as.integer(task_id), isTest=test)
x
