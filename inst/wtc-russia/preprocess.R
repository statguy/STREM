#library(devtools); install_github("statguy/Winter-Track-Counts")
library(parallel)
library(doMC)
registerDoMC(cores=detectCores())
library(WTC)
source("~/git/Winter-Track-Counts/setup/WTC-Boot.R")

context <- Context$new(resultDataDirectory=wd.data.results, processedDataDirectory=wd.data.processed, rawDataDirectory=wd.data.raw, scratchDirectory=wd.scratch, figuresDirectory=wd.figures, verbose=TRUE)
study <- RussiaWTCStudy$new(context=context, response="canis.lupus")
study$preprocess()
studyArea <- RussiaStudyArea$new(context=context)$saveBoundary(study=study)
study <- FinlandRussiaWTCStudy$new(context=context, response="canis.lupus")
study$preprocess()

if (F) {
  intersections <- study$loadIntersections()
  study$studyArea$habitat <- raster(extent(c(2.2,7.5,7.4,11.2)*1e6), nrows=200, ncols=600, crs=study$studyArea$proj4string)
  values(study$studyArea$habitat) <- 0
  plot(study$studyArea$boundary)
  plot(extent(study$studyArea$habitat), add=T)
  points(coordinates(intersections$intersections))
  points(coordinates(subset(intersections$intersections, district=="Finland")), col="red")
}