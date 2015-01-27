# Run test:
# ./parallel_r.py -t 1:50 -n 60 -l 10.0 -b ~/tmp/blacklist.txt -v ~/git/Winter-Track-Counts/inst/simulation/population_size.R test A SmoothModel-nbinomial-matern-ar1
# Run full:
# ./parallel_r.py -t 1:50 -n 60 -l 10.0 -b ~/tmp/blacklist.txt -v ~/git/Winter-Track-Counts/inst/simulation/population_size.R notest A SmoothModel-nbinomial-matern-ar1
# ./parallel_r.py -t 1:50 -n 60 -l 10.0 -b ~/tmp/blacklist.txt -v ~/git/Winter-Track-Counts/inst/simulation/population_size.R notest A SmoothModel-nbinomial-ar1
# ./parallel_r.py -t 1:50 -n 60 -l 10.0 -b ~/tmp/blacklist.txt -v ~/git/Winter-Track-Counts/inst/simulation/population_size.R notest A FMPModel
# ./parallel_r.py -t 1:5 -n 6 -l 10.0 -b ~/tmp/blacklist.txt -v ~/git/Winter-Track-Counts/inst/simulation/population_size.R notest Acombined SmoothModel-nbinomial-matern-ar1
# ./parallel_r.py -t 1:5 -n 6 -l 10.0 -b ~/tmp/blacklist.txt -v ~/git/Winter-Track-Counts/inst/simulation/population_size.R notest Acombined SmoothModel-nbinomial-ar1
# ./parallel_r.py -t 1:5 -n 6 -l 10.0 -b ~/tmp/blacklist.txt -v ~/git/Winter-Track-Counts/inst/simulation/population_size.R notest Acombined FMPModel
#
# R --vanilla --args notest E SmoothModel-nbinomial-ar1 3 < ~/git/Winter-Track-Counts/inst/simulation/population_size.R

# library(devtools); install_github("statguy/Winter-Track-Counts")
# echo 'library(devtools); install_github("statguy/Winter-Track-Counts")' | R --slave
# args <- c("notest","A","1")

if (F) {
scenario<-"Ecombined"
isTest<-F
modelName<-"SmoothModel-nbinomial-ar1"
iteration<-as.integer(1)
}

population_size <- function(scenario, modelName, iteration, isTest, otherTest=F) {
  readHabitatIntoMemory <- if (substr(scenario, 1, 1) == "E" || substr(scenario, 1, 1) == "F") TRUE else FALSE
  mss <- getMSS(scenario=scenario, isTest=isTest, readHabitatIntoMemory=readHabitatIntoMemory)
  study <- mss$study
  
  if (otherTest) {
    #iteration <- as.integer(50)
    estimates <- study$getModel(modelName=modelName, iteration=iteration)
    estimates <- study$loadEstimates(estimates=estimates)
    estimates$collectEstimates()
    
    #estimates$collectHyperparameters()
    populationSize <- estimates$getPopulationSize(withHabitatWeights=mss$hasHabitatWeights())
    populationSize$loadValidationData()
    populationSize
    colSums(populationSize$sizeData[,-1])
    colMeans(populationSize$sizeData[,-1])
    populationSize$plotPopulationSize()
    
    #study$loadPopulationSize(iteration=iteration, modelName="SmoothModel-nbinomial-ar1")
  }
  else {
    estimates <- study$getModel(modelName=modelName, iteration=iteration)
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
modelName <- extraArgs[1]
population_size(scenario=scenario, modelName=modelName, iteration=as.integer(task_id), isTest=isTest)
