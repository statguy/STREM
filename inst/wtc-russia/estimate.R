#library(devtools); install_github("statguy/Winter-Track-Counts")
library(parallel)
library(doMC)
registerDoMC(cores=detectCores())
library(WTC)
source("~/git/Winter-Track-Counts/setup/WTC-Boot.R")

if (test) {
  context <- Context$new(resultDataDirectory=wd.data.results, processedDataDirectory=wd.data.processed, rawDataDirectory=wd.data.raw, scratchDirectory=wd.scratch, figuresDirectory=wd.figures, verbose=TRUE)
  study <- FinlandRussiaWTCStudy$new(context=context, response="canis.lupus")
  
  intersections <- study$loadIntersections()
  intersections$intersections <- subset(intersections$intersections, year==2005)
  
  #meshParams <- list(maxEdge=c(1e6, 3e6), cutOff=0.5e6, coordsScale=1e-6)
  meshParams <- list(maxEdge=c(.1e6, 3e6), cutOff=0.15e6, coordsScale=1e-6)
  model <- FinlandRussiaSmoothModel$new(study=study)
  model$setup(intersections=intersections, meshParams=meshParams, offsetScale=1000^2)
  
  model$estimate(save=TRUE)
  
  populationDensity <- study$getPopulationDensity()
  plot(populationDensity$mean$rasterStack[[1]], col=terrain.colors(99))
}
else {
  estimate <- function(response) {
    context <- Context$new(resultDataDirectory=wd.data.results, processedDataDirectory=wd.data.processed, rawDataDirectory=wd.data.raw, scratchDirectory=wd.scratch, figuresDirectory=wd.figures, verbose=TRUE)
    study <- FinlandRussiaWTCStudy$new(context=context, response=response)
    
    meshParams <- list(maxEdge=c(.1e6, .3e6), cutOff=.08e6, coordsScale=1e-6)
    interceptPriorParams <- list(mean=200, sd=199)
    model <- study$estimate(meshParams=meshParams, interceptPriorParams=interceptPriorParams, save=TRUE)
  }
}
