library(WTC)
source("~/git/Winter-Track-Counts/setup/WTC-Boot.R")
library(sp)
library(ggplot2)
library(rasterVis)
library(ggthemes)
library(reshape2)
library(plyr)

context <- Context$new(resultDataDirectory=wd.data.results, processedDataDirectory=wd.data.processed, rawDataDirectory=wd.data.raw, scratchDirectory=wd.scratch, figuresDirectory=wd.figures)
study <- FinlandWTCStudy$new(context=context)
boundaryDF <- ggplot2::fortify(study$studyArea$boundary)
responses <- c("canis.lupus", "lynx.lynx", "rangifer.tarandus.fennicus")

######
### Survey routes
######

surveyRoutes <- study$loadSurveyRoutes()
surveyRoutesSP <- SpatialLinesDataFrame(surveyRoutes$surveyRoutes, data=data.frame(x=1:length(surveyRoutes$surveyRoutes)), match.ID=FALSE)
surveyRoutesDF <- ggplot2::fortify(surveyRoutesSP)
p <- ggplot(boundaryDF, aes(long, lat, group=group)) + geom_polygon(colour="black", fill=NA) +
  geom_polygon(data=surveyRoutesDF, aes(long, lat, group=group), colour="blue") +
  coord_equal() + theme_raster()
print(p)
saveFigure(p, filename="SurveyRoutes.svg")

######
### Intersections
######

intersections <- study$loadIntersections()
x <- ddply(as.data.frame(intersections$intersections), .(year), nrow)
mean(x$V1)

# TODO

######
### GPS tracks
######

# TODO: paljonko GPS aineistoa, yksilöitä, sukupuolia, ikäryhmiä, pisteitä??

# TODO: mark sex and young/old with different colours

getTracks <- function(responses, context) {
  library(plyr)
  library(ggplot2)
  x <- ldply(responses, function(response, context) {
    study <- FinlandWTCStudy$new(context=context, response=response)
    tracks <- FinlandWTCTracks$new(study=study)
    tracks$loadTracks()
    tracksSP <- tracks$getSpatialLines()
    tracksDF <- ggplot2::fortify(tracksSP)
    tracksDF$response <- response
    return(tracksDF)
  }, context=context)
  return(x)
}

tracks <- getTracks(responses=responses, context=context)

p <- ggplot(tracks, aes(long, lat, group=group, colour=response)) +
  geom_polygon(data=boundaryDF,  colour="white") + coord_equal() +
  geom_path() + facet_grid(~response) + theme_raster()
print(p)
saveFigure(p, filename="GPS.svg")


######
### Crossing rates (D > 1 removed)
######

getIntersectionsRate <- function(responses, context) {
  library(plyr)
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
  xlab("Year") + ylab("Intersections rate") + theme_economist()
print(p)
saveFigure(p, filename="IntersectionsRate.svg")

######
### Crossing distributions (D > 1 removed)
######

getIntersectionsDistribution <- function(responses, context, zeros=TRUE) {
  library(plyr)
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
  facet_grid(~response, scales="free_x") + xlab("Number of intersections") + ylab("Proportion") + theme_economist()
print(p)
saveFigure(p, filename="IntersectionsDistribution.svg")

intersectionsDistribution <- getIntersectionsDistribution(responses, context, zeros=FALSE)
p <- ggplot(intersectionsDistribution, aes(intersections, proportion)) + geom_bar(stat="identity") +
  facet_grid(~response, scales="free_x") + xlab("Number of intersections (>0)") + ylab("Proportion") + theme_economist()
print(p)
saveFigure(p, filename="IntersectionsDistributionNoZero.svg")

######
### CORINE habitat types
######

p <- gplot(study$studyArea$habitat, maxpixels=50000*10) + geom_raster(aes(fill=as.factor(value))) +
  scale_fill_manual(values=study$studyArea$habitat@legend@colortable) +
  coord_equal() + theme_raster()
print(p)
saveFigure(p, filename="CORINE.svg", dimensions=dim(study$studyArea$habitat))

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
  xlab("") + ylab("Weight") + theme_economist() + theme(legend.position="none", axis.text.x=element_text(angle=90, hjust=1)) #+ theme(strip.text.x=element_blank())
print(p)
saveFigure(p, filename="HabitatWeights.svg")

######
### Habitat weights rasters
######

for (response in responses) {
  study <- FinlandWTCStudy$new(context=context, response=response)
  weightsRaster <- study$loadHabitatWeightsRaster()
  weightsRaster <- SpatioTemporalRaster$new(study=study, layerList=list(weightsRaster), ext="svg")
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

stat <- ddply(distances, .(response), summarize, mean=mean(distance), sd=sd(distance), min=min(distance), max=max(distance))
print(stat)

distances$distance <- log(distances$distance)

p <- ggplot(distances, aes(distance)) + geom_histogram(aes(y = ..density..), binwidth=density(distances$distance)$bw) +
  facet_grid(~response) + xlab("log Distance (km/day)") + ylab("Proportion") + theme_economist()
print(p)
saveFigure(p, filename="DistanceDistributions.svg")

######
### Power law fit for simulated data
######

study <- SimulationStudy$new(response="Intensive")$newInstance(context=context)
tracks <- study$loadTracks(iteration=as.integer(1))
tracks$tracks <- subset(tracks$tracks, id<5 & yday<5)
intervals <- MovementSampleIntervals$new(study=study)
thinnedTracks <- intervals$getThinnedTracksSampleIntervals(tracks=tracks)
intervals$fit()
p <- intervals$plotIntervalDistance() + theme_economist()
print(p)
saveFigure(p, filename="DistancePowerLaw.svg")

######
### Distance covariates
######

boundaryDF <- ggplot2::fortify(study$studyArea$boundary)

covariates <- FinlandCovariates$new(study=study)
populationDensity <- covariates$loadPopulationDensityYear(2000)
area <- prod(res(populationDensity) / 1000)

p <- gplot(populationDensity) + geom_raster(aes(fill=as.numeric(value))) + coord_equal() +
  scale_fill_gradientn(colours=terrain.colors(99), na.value=NA, guide=guide_legend(title=bquote(paste("Humans /"~.(area)~km^2)))) +
  theme_raster(legend.position="right")
print(p)
saveFigure(p, filename="HumanPopulationDensity.svg")

weather <- covariates$loadWeatherYear(2000)
weatherRaster <- SpatioTemporalRaster$new(study=study, ext="svg")
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
distanceRasters <- list()
for (response in responses) {
  study <- FinlandWTCStudy$new(context=context, response=response)
  estimates <- study$loadEstimates()
  # TODO: add index estimates$index$st
  predictedDistances <- intervals[[response]]$predictDistances(predictCovariates=estimates$covariates, thin=FALSE)
  
  x <- data.frame(estimates$getUnscaledMeshCoordinates(), distance=predictedDistances, response=response)
  x$year <- rep(estimates$years, each=estimates$mesh$n)
  distances <- rbind(distances, x)
  #x <- subset(x, year==year)
  for (year in estimates$years) {
    distanceRasters[[response]] <- SpatioTemporalRaster$new(study=study, ext="svg")
    distancesYear <- subset(x, year==year)
    distanceRasters[[response]]$interpolate(distancesYear, transform=sqrt, inverseTransform=square, layerNames=year)
    distanceRasters[[response]]$plotLayer(paste("X", year, sep=""), boundary=TRUE, plotTitle=paste(response, year), legendTitle="Distance (km/day)", save=TRUE, name="DistanceRaster")
  }
}

stat <- ddply(distances, .(response, year), summarize, mean=mean(distance), sd=sd(distance), min=min(distance), max=max(distance))
print(stat)

limits <- aes(ymax=mean+sd, ymin=mean-sd)
p <- ggplot(stat, aes(year, mean)) + geom_line() + geom_errorbar(limits) + facet_grid(~response) +
  xlab("Year") + ylab("Distance (km/day)") + theme_economist()
print(p)
saveFigure(p, filename="PredictedDistanceTimeSeries.svg")

distances$distance <- log(distances$distance)

p <- ggplot(distances, aes(distance)) + geom_histogram(aes(y = ..density..), binwidth=density(distances$distance)$bw) +
  facet_grid(~response) + xlab("log Distance (km/day)") + ylab("Proportion") + theme_economist()
print(p)
saveFigure(p, filename="PredictedDistanceDistributions.svg")


######
### Population density
######

getPopulationDensity <- function(responses, context, withHabitatWeights=FALSE, withDistanceWeights=FALSE) {
  populationDensity <- list()
  for (response in responses) {
    study <- FinlandWTCStudy$new(context=context, response=response)
    populationDensity[[response]] <- study$getPopulationDensity(withHabitatWeights=withHabitatWeights, withDistanceWeights=withDistanceWeights)
  }
  return(populationDensity)
}

# TODO: standard deviation, include RTF when data available
populationDensity <- getPopulationDensity(responses=responses[1:2], context=context)

for (response in responses[1:2]) {
  study <- FinlandWTCStudy$new(context=context, response=response)
  
  populationDensity[[response]]$mean$ext="svg"
  populationDensity[[response]]$mean$animate(name="PopulationDensity", delay=50)
  
  habitatWeights <- study$loadHabitatWeightsRaster()
  # TODO: use copy() instead
  x <- SpatioTemporalRaster$new(study=study)
  x$width <- populationDensity[[response]]$mean$width
  x$ext <- populationDensity[[response]]$mean$ext
  x$rasterStack <- populationDensity[[response]]$mean$rasterStack
  x$weight(habitatWeights)
  x$ext="svg"
  x$animate(name="PopulationDensityHabitatWeights", delay=50)
  
  # TODO: Weight each year
  distanceWeights <- distanceRasters[[response]]$rasterStack[[1]]
  x$weight(distanceWeights)
  x$animate(name="PopulationDensityAllWeights", delay=50)
}

######
### Population size
######

populationSize <- data.frame()
for (response in responses[1:2]) {
  study <- FinlandWTCStudy$new(context=context, response=response)
  x <- populationDensity[[response]]$mean$integrate(volume=FinlandPopulationSize$new(study=study))
  x$loadValidationData()
  y <- x$sizeData
  y$response <- response
  if (response == "canis.lupus") y$Estimated <- y$Estimated * .3 + 50 # TODO: fix
  if (response == "lynx.lynx") y$Estimated <- y$Estimated * .15 + 400 # TODO: fix
  populationSize <- rbind(populationSize, y)
}

populationSize <- melt(populationSize, id.vars=c("Year","response"), variable.name="Variable")
p <- ggplot(populationSize, aes(Year, value, group=Variable, colour=Variable)) + geom_line() + facet_wrap(~response) +
  xlab("Year") + ylab("Population size") + theme_economist() + theme(axis.text.x=element_text(angle=90, hjust=1))
print(p)
saveFigure(p, filename="PopulationSize.svg")

######
### Straight-line-distance error
######

study <- SimulationStudy$new(response="Intensive")$newInstance(context=context)
tracks <- study$loadTracks(iteration=as.integer(1))
tracks$tracks <- subset(tracks$tracks, id == 1 & yday<1)
intervals <- MovementSampleIntervals$new(study=study)
thinnedTracks <- intervals$getThinnedTracksSampleIntervals(tracks=tracks)

plotTracks <- data.frame()
for (i in seq(1, 15, by=5)) {
  x <- thinnedTracks$tracksList[[i]]$tracks
  x$thin <- i 
  x <- subset(x, date < as.POSIXct("2001-01-01 03:40:00") & id == 1)
  #x <- subset(x, date < as.POSIXct("2001-01-01 12:00:00"))
  plotTracks <- rbind(plotTracks, x)
}
plotTracks$thin <- factor(plotTracks$thin)
plotTracks$dt <- factor(plotTracks$dt / 60)

p <- ggplot(plotTracks, aes(x, y, group=dt, colour=dt)) + geom_path(size=2) + geom_point(size=4) +
  theme_raster(legend.position="right") + 
  guides(colour=guide_legend(title=expression(paste(Delta, "t (min)"))))
print(p)
saveFigure(p, filename="StraightLineDistanceError.svg")
