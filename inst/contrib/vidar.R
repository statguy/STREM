library(parallel)
library(doMC)
registerDoMC(cores=detectCores())
library(WTC)
source("~/git/Winter-Track-Counts/setup/WTC-Boot.R")

response <- "lynx.lynx"
context <- Context$new(resultDataDirectory=wd.data.results, processedDataDirectory=wd.data.processed, rawDataDirectory=wd.data.raw, scratchDirectory=wd.scratch, figuresDirectory=wd.figures)
study <- FinlandWTCStudy$new(context=context, response=response, distanceCovariatesModel=~populationDensity+rrday+snow+tday-1, trackSampleInterval=2)
model <- FinlandSmoothModelSpatioTemporal(study=study)$setModelName("nbinomial", paste("matern", "ar1", sep="-"))
model$offsetScale <- 1000^2 # quickfix
model$loadEstimates()
model$collectEstimates()

library(raster)

grid <- raster(extent(320 * 1e4, 375 * 1e4, 665 * 1e4, 775 * 1e4), nrows=22, ncols=11)
grid[] <- rep(c(0,.5), times=ncell(grid)/2)
plot(grid)
plot(study$studyArea$boundary, add=T)
x <- data.frame(model$getUnscaledObservationCoordinates(), mean=model$data$fittedMean, year=model$data$year)

# Interpolate density from estimates
density <- study$getPopulationDensity(model=model, withHabitatWeights=TRUE)

# Save original density
densityOriginal <- writeRaster(density$mean$rasterStack, filename="~/tmp/lynx_density", overwrite=T)

# Resample density
densityResampled <- resample(density$mean$rasterStack, grid, filename="~/tmp/lynx_density_resampled", overwrite=T)

# Center value density
# TODO: handle NA's
y <- grid
x <- density$mean$rasterStack
ag <- res(y)/res(x)
densityCenterValue <- aggregate(x, ag, function(x, na.rm) x[length(x)/2], filename="~/tmp/lynx_density_centervalue", overwrite=T)

r1 <- raster("~/tmp/lynx_density")
plot(r1[[10]])
r2 <- raster("~/tmp/lynx_density_resampled")
plot(r2[[10]])
r3 <- raster("~/tmp/lynx_density_centervalue")
plot(r3[[10]])
