# Run test:
# ./parallel_r.py -t 1:50 -n 60 -l 10.0 -b ~/tmp/blacklist.txt -v ~/git/Winter-Track-Counts/inst/simulation/population_size.R test A
# Run full:
# ./parallel_r.py -t 1:50 -n 60 -l 10.0 -b ~/tmp/blacklist.txt -v ~/git/Winter-Track-Counts/inst/simulation/population_size.R notest A
# R --vanilla --args notest A 1 < population_size.R

# library(devtools); install_github("statguy/Winter-Track-Counts")
# echo 'library(devtools); install_github("statguy/Winter-Track-Counts")' | R --slave
# args <- c("notest","A","1")


population_size <- function(scenario, iteration, isTest, otherTest=F) {
  mss <- getMSS(scenario=scenario, isTest=isTest)
  study <- mss$study
  
  if (otherTest) {
    #iteration <- as.integer(50)
    estimates <- SimulatedSmoothModelSpatioTemporal(study=study, iteration=iteration)
    estimates <- study$loadEstimates(estimates=estimates)
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

library(parallel)
library(doMC)
registerDoMC(cores=detectCores())
library(WTC)
source("~/git/Winter-Track-Counts/setup/WTC-Boot.R")

parseArguments()
x <- population_size(scenario=scenario, iteration=as.integer(task_id), isTest=isTest)
x
