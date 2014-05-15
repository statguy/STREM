library(parallel)
library(doMC)
registerDoMC(cores=detectCores())
library(WTC)
source("~/git/Winter-Track-Counts/setup/WTC-Boot.R")

context <- Context$new(resultDataDirectory=wd.data.results, processedDataDirectory=wd.data.processed, rawDataDirectory=wd.data.raw, scratchDirectory=wd.scratch, figuresDirectory=wd.figures)
mss <- if (scenario == "A") MovementSimulationScenarioA$new()$newInstance(context=context) else stop("Unknown scenario ", scenario)
study <- mss$study

iterations <- study$context$getIterationIds(dir=study$context$resultDataDirectory, name="PopulationSize", response=study$response, region=study$studyArea$region)
if (length(iterations) == 0)
  stop("No files found.")

library(plyr)
populationSize <- ldply(iterations, function(iteration, study) {
  x <- study$loadPopulationSize(iteration=iteration)
  return(data.frame(x$sizeData, iteration=iteration))
}, study=study)
populationSize <- arrange(populationSize, iteration, Year)
print(populationSize)
#populationSize <- subset(populationSize, Estimated < 10000) # Some results strange... fixed with prior

ddply(populationSize, .(Year), summarise, ObservedMean=mean(Observed), ObservedSD=sd(Observed), EstimatedMean=mean(Estimated), EstimatedSD=sd(Estimated))
ddply(populationSize, .(iteration), summarise, Observed=mean(Observed, na.rm=T), Estimated=mean(Estimated, na.rm=T), iteration=mean(iteration))

library(ggplot2)
library(reshape2)

x <- melt(populationSize, id.vars="Year")
p <- ggplot(x, aes(x=Year, y=value, group=variable, colour=variable)) +
  geom_point(stat="summary", stat_params=list(fun.y="mean")) +
  geom_line(stat="summary", stat_params=list(fun.y="mean")) +
  stat_summary(fun.data=mean_se, geom="errorbar") +
  ylab("Population size")
print(p)
ggsave(study$context$getFileName(dir=study$context$figuresDirectory, name="PopulationSize", response=study$response, region=study$studyArea$region, ext=".png"), p, width=10, height=8)

###

intersectionsDistribution <- llply(iterations, function(iteration, study) {
  x <- study$loadIntersections(iteration=iteration)
  return(x$intersections$intersections)
}, study=study)
round(prop.table(table(Reduce(c, intersectionsDistribution)))*100,2)
