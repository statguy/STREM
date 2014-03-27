library(WTC)
source("~/git/Winter-Track-Counts/setup/WTC-Boot.R")
library(ggplot2)
library(rasterVis)

context <- Context$new(resultDataDirectory=wd.data.results, processedDataDirectory=wd.data.processed, rawDataDirectory=wd.data.raw, scratchDirectory=wd.scratch, figuresDirectory=wd.figures)
study <- FinlandWTCStudy$new(context=context)
responses <- c("canis.lupus", "lynx.lynx", "rangifer.tarandus.fennicus")

######
### GPS tracks
######

# TODO: mark sex and young/old with different colours

getTracks <- function(responses, context) {
  library(plyr)
  library(ggplot2)
  x <- ldply(responses, function(response, context) {
    study <- FinlandWTCStudy$new(context=context, response=response)
    tracks <- FinlandWTCTracks$new(study=study)
    tracks$loadTracks()
    tracksSP <- tracks$getSpatialLines()
    tracksDF <- fortify(tracksSP)
    tracksDF$response <- response
    return(tracksDF)
  }, context=context)
  return(x)
}

boundaryDF <- fortify(study$studyArea$boundary)
tracks <- getTracks(responses=responses, context=context)

p <- ggplot(tracks, aes(long, lat, group=group, colour=response)) +
  geom_polygon(data=boundaryDF,  colour="white") + coord_equal() +
  geom_path() + facet_grid(~response) + theme_raster()
print(p)
ggsave(p, filename=file.path(context$figuresDirectory, "GPS.pdf"), width=8)

######
### Crossing rates (D > 1 removed)
######

getIntersectionsRate <- function(responses, context) {
  x <- ldply(responses, function(response, context) {
    study <- FinlandWTCStudy$new(context=context, response=response)
    intersections <- study$loadIntersections()
    x <- ddply(intersections$intersections@data, .(year), function(x, response) {
      data.frame(intersections=sum(x$intersections / (x$length/1000) / x$duration) / nrow(x), response=response)
    }, response=response)
    return(x)
  }, context=context)
  return(x)
}

intersectionsRate <- getIntersectionsRate(responses, context)
p <- ggplot(intersectionsRate, aes(year, intersections)) + geom_line() + facet_grid(~response) +
  xlab("Year") + ylab("Intersections rate") + theme_bw()
print(p)
ggsave(p, filename=file.path(context$figuresDirectory, "IntersectionsRate.pdf"), width=8)


######
### Crossing distributions (D > 1 removed)
######

getIntersectionsDistribution <- function(responses, context, zeros=TRUE) {
  x <- ldply(responses, function(response, context, zeros) {
    study <- FinlandWTCStudy$new(context=context, response=response)
    intersections <- study$loadIntersections()
    intersections <- if (zeros == TRUE) table(intersections$intersections@data$intersections)
    else table(intersections$intersections@data$intersections[intersections$intersections@data$intersections>0])
    x <- as.data.frame(intersections)
    x$Freq <- x$Freq / sum(x$Freq)
    names(x) <- c("intersections", "proportion")
    x$intersections <- as.integer(as.character(x$intersections))
    x$response <- response
    return(x)
  }, context=context, zeros=zeros)
  return(x)
}

intersectionsDistribution <- getIntersectionsDistribution(responses, context)
p <- ggplot(intersectionsDistribution, aes(intersections, proportion)) + geom_bar(stat="identity") +
  facet_grid(~response, scales="free_x") + xlab("Number of intersections") + ylab("Proportion") + theme_bw()
print(p)
ggsave(p, filename=file.path(context$figuresDirectory, "IntersectionsDistribution.pdf"), width=8)

intersectionsDistribution <- getIntersectionsDistribution(responses, context, zeros=FALSE)
p <- ggplot(intersectionsDistribution, aes(intersections, proportion)) + geom_bar(stat="identity") +
  facet_grid(~response, scales="free_x") + xlab("Number of intersections (>0)") + ylab("Proportion") + theme_bw()
print(p)
ggsave(p, filename=file.path(context$figuresDirectory, "IntersectionsDistributionNoZero.pdf"), width=8)

######
### CORINE habitat types
######

p <- gplot(study$studyArea$habitat, maxpixels=50000*10) + geom_raster(aes(fill=as.factor(value))) +
  scale_fill_manual(values=study$studyArea$habitat@legend@colortable) +
  coord_equal() + theme_raster()
print(p)
aspectRatio <- dim(study$studyArea$habitat)[2] / dim(study$studyArea$habitat)[1]
ggsave(p, filename=file.path(context$figuresDirectory, "CORINE.pdf"), width=8, height=8 * aspectRatio)

######
### Habitat weights
######

getHabitatWeights <- function(responses, context) {
  library(plyr)
  x <- ldply(responses, function(response, context) {
    study <- FinlandWTCStudy$new(context=context, response=response)
    weights <- study$loadHabitatWeights()$getHabitatSelectionWeights()
    n <- names(weights)
    print(n)
    weights <- data.frame(habitat=factor(n, levels=n), weights)
    weights$response <- response
    return(weights)
  }, context=context)
  return(x)
}

weights <- getHabitatWeights(responses=responses, context=context)
p <- ggplot(weights, aes(habitat, weights, fill=habitat)) +
  geom_bar(stat="identity") + facet_grid(~response) + scale_fill_manual(values=c("#beaed4","#ffff99","#7fc97f","#fdc086","#386cb0")) +
  xlab("") + ylab("Weight") + theme_bw()
print(p)
ggsave(p, filename=file.path(context$figuresDirectory, "HabitatWeights.pdf"), width=8)

######
### Habitat weights rasters
######

for (response in responses) {
  study <- FinlandWTCStudy$new(context=context, response=response)
  weightsRaster <- study$loadHabitatWeightsRaster()
  weightsRaster <- SpatioTemporalRaster$new(study=study, layerList=list(weightsRaster), ext="pdf")
  weightsRaster$plotLayer(layerName=1, save=TRUE, name="HabitatWeights")
}

######
### Distance distributions
######

distances <- list()
for (response in responses) {
  study <- FinlandWTCStudy$new(context=context, response=response)
  tracks <- study$loadTracks()
  intervals <- tracks$getSampleIntervals()
  distances <- rbind(distances, data.frame(distance=intervals$intervals$dist, response=response))  
}

stat <- ddply(distances, .(response), function(x) {
  data.frame(mean=mean(x$distance), sd=sd(x$distance))
})
print(stat)

p <- ggplot(distances, aes(distance)) + geom_histogram(aes(y = ..density..), binwidth=density(distances$distance)$bw) +
  facet_grid(~response) + xlab("Distance (km/day)") + ylab("Proportion") + theme_bw()
print(p)
ggsave(p, filename=file.path(context$figuresDirectory, "DayDistances.pdf"), width=8)

######
### Power law fit for simulated data
######

study <- SimulationStudy$new(response="Intensive")$newInstance(context=context)
tracks <- study$loadTracks(iteration=as.integer(1))
tracks$tracks <- subset(tracks$tracks, id<5 & yday<5)
intervals <- MovementSampleIntervals$new(study=study)
thinnedTracks <- intervals$getThinnedTracksSampleIntervals(tracks=tracks)
intervals$fit()
p <- intervals$plotIntervalDistance() + theme_bw()
print(p)
ggsave(p, filename=file.path(context$figuresDirectory, "DistancePowerLaw.pdf"), width=8)

######
### Distance covariates
######

boundaryDF <- fortify(study$studyArea$boundary)

covariates <- FinlandCovariates$new(study=study)
populationDensity <- covariates$loadPopulationDensityYear(2000)
area <- prod(res(populationDensity) / 1000)

p <- gplot(populationDensity) + geom_raster(aes(fill=as.numeric(value))) + coord_equal() +
  scale_fill_gradientn(colours=terrain.colors(99), na.value=NA, guide=guide_legend(title=bquote(paste("Humans /"~.(area)~km^2)))) +
  theme_raster(legend.position="right")
print(p)
ggsave(p, filename=file.path(context$figuresDirectory, "HumanPopulationDensity.pdf"), width=8)


x <- covariates$loadWeatherYear(2000)
gplot(x$month1snow) + geom_raster(aes(fill=value)) + coord_equal() +
  scale_fill_gradientn(colours=terrain.colors(20), na.value=NA, guide=guide_legend(title="Snow cover")) +
  theme_raster(legend.position="right") +
  geom_polygon(data=boundaryDF, aes(long, lat), colour="white", fill=NA)

weather <- covariates$loadWeatherYear(2000)
weatherRaster <- SpatioTemporalRaster$new(study=study)
weatherRaster$addLayer(weather$month1snow, "snowcover")
weatherRaster$addLayer(weather$month1rrday, "precip")
weatherRaster$addLayer(weather$month1tday, "temp")
weatherRaster$plotLayer("snowcover", boundary=T, plotTitle="January 2000", legendTitle="Snow depth (cm)", save=TRUE, name="Weather")
weatherRaster$plotLayer("precip", boundary=T, plotTitle="January 2000", legendTitle="Precipitation (?)", save=TRUE, name="Weather")
weatherRaster$plotLayer("temp", boundary=T, plotTitle="January 2000", legendTitle="Temperature (°C)", save=TRUE, name="Weather")

######
### Predicted distance histograms and rasters
######

intervals <- list()
for (response in responses) {
  study <- FinlandWTCStudy$new(context=context, response=response)
  intervals[[response]] <- study$getSampleIntervals()
  tracks <- study$loadTracks()
  intervals[[response]]$getThinnedTracksSampleIntervals(tracks=tracks)
}

distances <- data.frame()
distanceRasters <- SpatioTemporalRaster$new(study=study, ext="pdf")
for (response in responses) {
  study <- FinlandWTCStudy$new(context=context, response=response)
  model <- log(distanceKm) ~ sqrt(sqrt(populationDensity)) + sqrt(rrday) + sqrt(snow) + log(intervalH) # TODO: problem with tday??
  intervals[[response]]$estimationResult <- glm(model, data=intervals[[response]]$getDailyDistanceData())
  print(intervals[[response]]$estimationResult)
  estimates <- study$loadEstimates()
  covariates <- estimates$covariates
  covariates$intervalH <- 2
  x <- exp(intervals[[response]]$predict(predictionData=covariates))
  x[x>50] <- 50
  x <- data.frame(estimates$getUnscaledMeshCoordinates(), distance=x, response=response)
  distances <- rbind(distances, x)
  distanceRasters$interpolate(x, transform=sqrt, inverseTransform=square, layerNames=response)
  distanceRasters$plotLayer(response, boundary=T, plotTitle=response, legendTitle="Distance (km/day)", save=TRUE, name="DistanceRaster")
}

stat <- ddply(distances, .(response), summarize, mean=mean(distance), sd=sd(distance), min=min(distance), max=max(distance))
print(stat)

p <- ggplot(distances, aes(distance)) + geom_histogram(aes(y = ..density..), binwidth=density(distances$distance)$bw) +
  facet_grid(~response) + xlab("Distance (km/day)") + ylab("Proportion") + theme_bw()
print(p)
ggsave(p, filename=file.path(context$figuresDirectory, "DistanceDistribution.pdf"), width=8)

######
### Straight-line-distance error
######

plotTracks <- data.frame()
for (i in seq(1, 15, by=5)) {
  x <- thinnedTracks$tracksList[[i]]$tracks
  x$thin <- i 
  x <- subset(x, date < as.POSIXct("2001-01-01 03:40:00"))
  plotTracks <- rbind(plotTracks, x)
}
plotTracks$thin <- factor(plotTracks$thin)
plotTracks$dt <- factor(plotTracks$dt / 60)

p <- ggplot(plotTracks, aes(x, y, group=dt, colour=dt)) + geom_path(size=2) + geom_point(size=4) +
  theme_raster(legend.position="right") + 
  guides(colour=guide_legend(title=expression(paste(Delta, "t (min)"))))
print(p)
ggsave(p, filename=file.path(context$figuresDirectory, "ThinnedTracks.pdf"), width=8, height=4)



### Unfinished: ###

######
### Population density and size
######

getPopulationDensity <- function(responses, context, withHabitatWeights=FALSE, withDistanceWeights=FALSE) {
  populationDensity <- list()
  for (response in responses) {
    study <- FinlandWTCStudy$new(context=context, response=response)
    populationDensity[[response]] <- study$getPopulationDensity(withHabitatWeights=withHabitatWeights, withDistanceWeights=withDistanceWeights)
  }
  return(populationDensity)
}

populationDensity <- getPopulationDensity(responses=responses, context=context)

for (response in responses[1:2]) {
  study <- FinlandWTCStudy$new(context=context, response=response)
  populationDensity[[response]]$mean$animate(name="PopulationDensity", delay=50)
}

for (response in responses[1:2]) {
  study <- FinlandWTCStudy$new(context=context, response=response)
  populationSize <- populationDensity[[response]]$mean$integrate(volume=FinlandPopulationSize$new(study=study))
  populationSize$loadValidationData()
}
