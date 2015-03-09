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
library(parallel)
library(doMC)
registerDoMC(cores=detectCores())
library(WTC)
source("~/git/Winter-Track-Counts/setup/WTC-Boot.R")

parseArguments()
  
scenario <- "E"
#scenario <- "A"
#scenario<-"Ecombined"
isTest<-F
#modelName<-"FMPModel"
#modelName<-"SmoothModelMean-nbinomial-ar1"
modelName<-"SmoothModelMean-nbinomial-ar1-priors1"
#modelName<-"SmoothModel-nbinomial-matern-ar1"
iteration<-as.integer(3)
readHabitatIntoMemory <- F
mss <- getMSS(scenario=scenario, isTest=isTest, readHabitatIntoMemory=FALSE)
study <- mss$study
#surveyRoutes <- mss$getSurveyRoutes()
#populationSize <- study$getPopulationSize2(modelName=modelName, iteration=iteration, save=F)
}

population_size <- function(scenario, modelName, iteration, isTest, otherTest=F) {
  readHabitatIntoMemory <- if (substr(scenario, 1, 1) == "E" || substr(scenario, 1, 1) == "F") TRUE else FALSE
  #mss <- getMSS(scenario=scenario, isTest=isTest, readHabitatIntoMemory=readHabitatIntoMemory)
  mss <- getMSS(scenario=scenario, isTest=isTest, readHabitatIntoMemory=FALSE)
  study <- mss$study
  surveyRoutes <- mss$getSurveyRoutes()
  
  if (otherTest) {
    study <- getMSS(scenario="Acombined")$study
    iteration <- as.integer(1)
    
    estimates <- study$getModel(modelName=modelName, iteration=iteration)
    estimates <- study$loadEstimates(estimates=estimates)
    estimates$collectEstimates()
    
    #estimates$collectHyperparameters()
    populationSize <- estimates$getPopulationSize()
    populationSize$loadValidationData()
    populationSize
    colSums(populationSize$sizeData[,-1])
    colMeans(populationSize$sizeData[,-1])
    populationSize$plotPopulationSize()
    
    #study$loadPopulationSize(iteration=iteration, modelName=modelName)
  }
  else {
    #iteration<-as.integer(6)
    
    populationSize <- study$getPopulationSize2(modelName=modelName, iteration=iteration, save=T)
    return(invisible(populationSize))
    
    #populationSize <- study$getPopulationSize(estimates=estimates, habitatWeights=habitatWeights)
    #return(invisible(populationSize))
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
