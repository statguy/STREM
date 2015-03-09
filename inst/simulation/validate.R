# ./parallel_r.py -t 1:50 -n 70 -l 10.0 -b ~/tmp/blacklist.txt -v ~/git/Winter-Track-Counts/inst/simulation/validate.R notest A SmoothModel-nbinomial-ar1
# ./parallel_r.py -t 1:50 -n 70 -l 10.0 -b ~/tmp/blacklist.txt -v ~/git/Winter-Track-Counts/inst/simulation/validate.R notest A SmoothModel-nbinomial-matern-ar1

# library(devtools); install_github("statguy/Winter-Track-Counts")

library(parallel)
library(doMC)
registerDoMC(cores=detectCores())
library(WTC)
source("~/git/Winter-Track-Counts/setup/WTC-Boot.R")

parseArguments()
modelName <- extraArgs[1]

nSamples <- 50
mss <- getMSS(scenario=scenario, isTest=isTest)
study <- mss$study
cutoff <- c(1,2000)
if (substr(scenario, 2, nchar(scenario)) == "combined") cutoff <- 10 * cutoff
validation <- Validation(study=study, populationSizeCutoff=cutoff)
iteration <- as.integer(task_id)

if (isTest) {
  populationSize <- validation$validateTemporalPopulationSize(modelName=modelName)
  populationSize
  print(validation$summarizePopulationSize(populationSize))
  print(summary(lm(Estimated~Observed, populationSize)))
  
  library(plyr)
  #iteration <- as.integer(1)
  iterations <- validation$getEstimatesFileIterations(modelName=modelName)
  populationSizeCI <- ldply(iterations, function(iteration) {
    validation$validateCredibilityIntervals(modelName=modelName, iteration=iteration, nSamples=nSamples, save=F)
  }, .parallel=T)
  print(validation$summarizePopulationSizeCI(populationSizeCI))
  print(validation$summarizePopulationSizeCI(populationSizeCI, probs=c(.25,.75)))
  
  
  x <- ddply(populationSizeCI, .(scenario, Year, iteration), function(x, probs) {
    y <- data.frame(Estimated=mean(x$Estimated), Observed=mean(x$Observed))
    q <- quantile(x$Estimated, probs=probs)
    y$Estimated.q1 <- q[1]
    y$Estimated.q2 <- q[2]
    return(y)
  }, probs=probs)
  xy <- merge(populationSize[,c("iteration","Estimated")], x, by="iteration")
  summary(lm(Estimated.y~Estimated.x, xy)) # OK !
  
} else {
  populationSizeCI <- validation$validateCredibilityIntervals(modelName=modelName, iteration=iteration, nSamples=nSamples, save=T)
}
