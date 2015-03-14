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

# R --vanilla --args notest E SmoothModel-nbinomial-matern-ar1 3 < ~/git/Winter-Track-Counts/inst/simulation/population_size.R
# R --vanilla --args notest E SmoothModel-nbinomial-matern-ar1 3 < ~/git/Winter-Track-Counts/inst/simulation/validation.R

# library(devtools); install_github("statguy/Winter-Track-Counts")
# echo 'library(devtools); install_github("statguy/Winter-Track-Counts")' | R --slave
# args <- c("notest","A","1")

if (F) {
library(parallel)
library(doMC)
registerDoMC(cores=round(detectCores()) * 3/4)
library(WTC)
source("~/git/Winter-Track-Counts/setup/WTC-Boot.R")

parseArguments()

scenario <- "E"
#scenario <- "A"
#scenario<-"Ecombined"
isTest<-F
#modelName<-"FMPModel"
#modelName<-"SmoothModelMean-nbinomial-ar1"
#odelName<-"SmoothModelMean-nbinomial-ar1-priors1"
#modelName<-"SmoothModel-nbinomial-matern-ar1"
iteration<-as.integer(1)
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
  surveyRoutes <- study$surveyRoutes
  
  if (otherTest) {
    library(WTC)
    source("~/git/Winter-Track-Counts/setup/WTC-Boot.R")
    scenario <- "F"
    readHabitatIntoMemory <- F
    modelName <- "SmoothModelMean-nbinomial-ar1-priors1"
    #modelName<-"SmoothModel-nbinomial-matern-ar1"
    study <- getMSS(scenario=scenario, readHabitatIntoMemory=FALSE)$study
    study$grassLocalTempDir <- grassLocalTempDir
    iteration <- as.integer(18)
    
    habitatWeights <- study$getHabitatWeights(iteration=iteration, readHabitatIntoMemory=FALSE)
    
    estimates <- study$getModel(modelName = modelName, iteration = iteration)
    estimates$loadEstimates()
    estimates$collectEstimates()    
    populationDensity <- estimates$getPopulationDensity(habitatWeights=habitatWeights, .parallel=FALSE)
    #habitatWeightsRaster <- habitatWeights$getWeightsRaster(grassLocalTempDir=study$grassLocalTempDir, save=TRUE)
    #writeRaster(habitatWeightsRaster, "public_html/x.tif", format="GTiff", overwrite=T)
    habitatWeightsRaster <- habitatWeights$getWeightsRaster()
    populationDensity$weight(habitatWeightsRaster)
    plot(populationDensity$rasterStack[[1]])     
    
    populationSize <- estimates$getPopulationSize(populationDensity, habitatWeightsRaster=habitatWeightsRaster)
    
    
  }
  else {
    study$grassLocalTempDir <- grassLocalTempDir
    populationSize <- study$getPopulationSize2(modelName=modelName, iteration=iteration, save=T)
    return(invisible(populationSize))
    
    #populationSize <- study$getPopulationSize(estimates=estimates, habitatWeights=habitatWeights)
    #return(invisible(populationSize))
  }
}


library(parallel)
library(doMC)
registerDoMC(cores=detectCores() * 1/2)
library(WTC)
source("~/git/Winter-Track-Counts/setup/WTC-Boot.R")

parseArguments()
modelName <- extraArgs[1]
population_size(scenario=scenario, modelName=modelName, iteration=as.integer(task_id), isTest=isTest)
