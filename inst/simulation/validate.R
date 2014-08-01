# library(devtools); install_github("statguy/Winter-Track-Counts")

library(parallel)
library(doMC)
registerDoMC(cores=detectCores())
library(WTC)
source("~/git/Winter-Track-Counts/setup/WTC-Boot.R")

parseArguments()

if (isTest) {
  nSamples <- 50
  populationSizeOverEstimate <- 200
} else {
  nSamples <- 50
  populationSizeOverEstimate <- 2000
}

modelName <- paste("SmoothModel", "nbinomial", "matern", "ar1", sep = "-")
mss <- getMSS(scenario=scenario, isTest=isTest)
study <- mss$study
validation <- Validation(study=study, populationSizeOverEstimate=populationSizeOverEstimate)

if (isTest) {
  populationSize <- validation$validateTemporalPopulationSize(modelName=modelName)
  populationSize
  print(validation$populationSizeSummary(populationSize))
  print(summary(lm(Estimated~Observed, populationSize)))
  
  #iteration <- as.integer(1)
  populationSizeCI <- validation$validateCredibilityIntervals(modelName=modelName, iteration=iteration, nSamples=nSamples, save=T)
  populationSizeCI
  
  ddply(populationSizeCI, .(scenario, Year), function(x, probs) {
    y <- data.frame(Estimated=mean(x$Estimated), Observed=mean(x$Observed))
    q <- quantile(x$Estimated, probs=probs)
    y$Estimated.q1 <- q[1]
    y$Estimated.q2 <- q[2]
    return(y)
  }, probs=c(.0025, .975))

} else {
  populationSize <- validation$validateTemporalPopulationSize(modelName=modelName)
  validation$populationSizeSummary(populationSize)
  summary(lm(Estimated~Observed, populationSize))  
  
  spatialCorrelation <- validation$validateSpatialPopulationSize(modelName=modelName)
  summary(lm(Correlation~True, spatialCorrelation))
  
  populationSizeCI <- validation$validateCredibilityIntervals(modelName=modelName, iteration=iteration, nSamples=nSamples, save=T)
  validation$loadCredibilityIntervalsValidation(modelName=modelName, iteration=iteration)
  
  ddply(populationSizeCI, .(scenario, Year), function(x, probs) {
    y <- data.frame(Estimated=mean(x$Estimated), Observed=mean(x$Observed))
    q <- quantile(x$Estimated, probs=probs)
    y$Estimated.q1 <- q[1]
    y$Estimated.q2 <- q[2]
    return(y)
  }, probs=c(.0025, .975))  
}

if (F) {
  summary(lm(Estimated~Observed, populationSize))
  p <- ggplot(size, aes(Observed, Estimated)) + geom_abline(colour="darkgrey") +
    geom_point() + geom_smooth(method="lm", colour="black") +
    facet_grid(~scenario)
 
  summary(lm(Correlation~True, spatialCorrelation))
  ggplot(spatialCorrelation, aes(True, Correlation)) + geom_point() + geom_smooth(method="lm", colour="black") +
    facet_grid(~scenario)
}
