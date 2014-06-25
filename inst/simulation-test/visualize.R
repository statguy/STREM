library(WTC)
source("~/git/Winter-Track-Counts/setup/WTC-Boot.R")
library(ggplot2)
library(rasterVis)


######
### Simulations (test study area)
######

if (F) {
  library(doMC)
  registerDoMC(cores=detectCores())
}
task_id <- 0; source(file.path(path.package("WTC"), "simulation-test", "simulate.R"))

simulate(scenario="A", nIterations=as.integer(1), plot=TRUE)
simulate(scenario="B", nIterations=as.integer(1), plot=TRUE)
simulate(scenario="C", nIterations=as.integer(1), plot=TRUE)
simulate(scenario="D", nIterations=as.integer(1), plot=TRUE)
simulate(scenario="E", nIterations=as.integer(1), plot=TRUE)
simulate(scenario="F", nIterations=as.integer(1), plot=TRUE)

context <- Context$new(resultDataDirectory=wd.data.results, processedDataDirectory=wd.data.processed, rawDataDirectory=wd.data.raw, scratchDirectory=wd.scratch, figuresDirectory=wd.figures)
study <- SimulationStudy$new()$setup(context=context, isTest=T)
png(file.path(context$figuresDirectory, "SimulationStudyArea-E-test.png"), bg="transparent")
plot(study$studyArea$habitat)
dev.off()

p <- ggplot(study$studyArea$toGGDF(), aes(long, lat, group=group)) + geom_path(size=1) + coord_equal() + theme_raster()
plot(p)
saveFigure(p, filename="SimulationStudyArea-A-test.svg", bg="transparent")
