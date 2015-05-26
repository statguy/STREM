#library(devtools); install_github("statguy/Winter-Track-Counts")
library(parallel)
library(doMC)
registerDoMC(cores=round(detectCores()))
library(WTC)
source("~/git/Winter-Track-Counts/setup/WTC-Boot.R")
library(reshape2)
library(boot)

parseArguments()
modelName <- extraArgs[1]
iteration <- as.integer(task_id)
#modelName <- "FMPModel"
#scenario <- "A"
#iteration <- as.integer(1)
mss <- getMSS(scenario=scenario)
study <- mss$study
context <- study$context
validation <- Validation(study=study, populationSizeCutoff=Inf)

nSamples <- 100

model <- study$getModel(modelName=modelName, iteration=iteration)
model$loadEstimates()
model$collectEstimates()
fitted <- model$getDensityEstimates()

#t <- length(model$years)
#n <- nrow(fitted) / t
#fitted$id <- rep(1:n, times=t)
#fitted <- dcast(fitted, year ~ id, value.var="density")

getPopulationSize <- function(fitted, index, model, readHabitatIntoMemory) {
  populationSize <- study$getPopulationSize(model, index=index, readHabitatIntoMemory=readHabitatIntoMemory, loadValidationData=F, save=F)
  return(populationSize$sizeData$Estimated)
}

findCI <- function(fitted, iteration, model, modelName) {
  populationSize <- study$loadPopulationSize(iteration, modelName)
  
  b <- ddply(fitted, .(year), function(x, iteration, model, populationSize) {
    message("Resampling year ", x$year[1], "...")
    year <- as.integer(x$year[1])
    true <- subset(populationSize$sizeData, Year == year)$Observed
    b <- boot(data=x$density, statistic=getPopulationSize, R=1000, model=model, readHabitatIntoMemory=F, parallel="multicore")
    #ci95 <- quantile(b$t, c(.025,.975))
    #ci50 <- quantile(b$t, c(.25,.75))
    bci <- boot.ci(b, conf=c(.95, .50), ype=c("basic"))
    if (!is.null(bci)) {      
      ci95 <- bci$basic[1,4:5]
      ci50 <- bci$basic[2,4:5]
      p95 <- ci95[1] <= true & ci95[2] >= true
      p50 <- ci50[1] <= true & ci50[2] >= true
      return(data.frame(year=year, p95=p95, p50=p50))
    }
    else return(data.frame(year=year, p95=NA, p50=NA))
  }, iteration=iteration, model=model, populationSize)
  return(b)
}

bootCI <- llply(1:nSamples, function(i) {
  message("Resampling ", i, "/", nSamples, "...")
  return(findCI(fitted, iteration, model, modelName))
}, .parallel=F)

bootCI <- do.call(rbind, bootCI)
save(bootCI, file=validation$getCredibilityIntervalsValidationFileName(modelName=modelName, iteration=iteration))
