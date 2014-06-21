library(WTC)
source("~/git/Winter-Track-Counts/setup/WTC-Boot.R")
library(sp)
library(ggplot2)
library(rasterVis)
library(ggthemes)
library(reshape2)
library(plyr)

context <- Context$new(resultDataDirectory=wd.data.results, processedDataDirectory=wd.data.processed, rawDataDirectory=wd.data.raw, scratchDirectory=wd.scratch, figuresDirectory=wd.figures)
study <- FinlandWTCStudy$new(context=context, response="canis.lupus")
study$studyArea$loadBoundary(thin=TRUE, tolerance=0.001)
boundaryDF <- study$studyArea$toGGDF()
responses <- c("canis.lupus", "lynx.lynx", "rangifer.tarandus.fennicus")

######
### Survey routes
######

surveyRoutes <- study$loadSurveyRoutes(findLengths=FALSE)
p <- ggplot(boundaryDF, aes(long, lat, group=group)) + geom_polygon(colour="black", fill=NA) +
  geom_polygon(data=surveyRoutes$toGGDF(), aes(long, lat, group=group), colour="blue", fill=NA) +
  coord_equal() + theme_raster()
print(p)
saveFigure(p, filename="SurveyRoutes.svg", bg="transparent")

p <- ggplot(boundaryDF, aes(long, lat, group=group)) + geom_polygon(colour="black", fill=NA, size=1) +
  geom_polygon(data=surveyRoutesDF, aes(long, lat, group=group), colour="blue", fill=NA, size=1) +
  coord_cartesian(c(3.605,3.638)*1e6, c(7.27,7.29)*1e6) + theme_raster(aspect.ratio=.6, panel.background=element_rect(fill="transparent", colour="gray", size=1))
print(p)
saveFigure(p, filename="SurveyRoutesZoom.svg", bg="transparent", width=3, height=3)

#surveyRouteDF <- ggplot2::fortify(surveyRoutesSP[1,])
p <- ggplot(subset(surveyRoutes$toGGDF(), id==1), aes(long, lat, group=group)) + geom_polygon(colour="blue", fill=NA) + coord_equal() + theme_raster()
print(p)
saveFigure(p, filename="SurveyRoute.svg", bg="transparent", width=2, height=2)

######
### Intersections
######

intersections <- study$loadIntersections(predictDistances=FALSE)
x <- ddply(as.data.frame(intersections$intersections), .(year), nrow)
mean(x$V1); sd(x$V1)

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
    tracksDF <- tracks$toGGDF(response)
    return(tracksDF)
  }, context=context)
  return(x)
}

tracks <- getTracks(responses=responses, context=context)

p <- ggplot(tracks, aes(long, lat, group=group, colour=response)) +
  geom_polygon(data=boundaryDF,  colour="white") + coord_equal() +
  geom_path() + facet_grid(~response) + theme_raster(20)
print(p)
saveFigure(p, filename="GPS.svg", bg="transparent")


######
### Crossing rates (D > 1 removed)
######

getIntersectionsRate <- function(responses, context) {
  library(plyr)
  x <- ldply(responses, function(response, context) {
    study <- FinlandWTCStudy$new(context=context, response=response)
    intersections <- study$loadIntersections(predictDistances=FALSE)
    x <- ddply(intersections$intersections@data, .(year), function(x, response) {
      data.frame(intersections=sum(x$intersections / (x$length/1000) / x$duration) / nrow(x), response=study$getPrettyResponse(response))
    }, response=response)
    return(x)
  }, context=context)
  return(x)
}

intersectionsRate <- getIntersectionsRate(responses, context)
p <- ggplot(intersectionsRate, aes(year, intersections)) + geom_line(size=1) + facet_grid(~response) +
  xlab("Year") + ylab("Intersections rate") + theme_presentation(16)
print(p)
saveFigure(p, filename="IntersectionsRate.svg", bg="transparent")

######
### Crossing distributions (D > 1 removed)
######

getIntersectionsDistribution <- function(responses, context, zeros=TRUE) {
  library(plyr)
  x <- ldply(responses, function(response, context, zeros) {
    study <- FinlandWTCStudy$new(context=context, response=response)
    intersections <- study$loadIntersections(predictDistances=FALSE)
    intersections <- if (zeros == TRUE) table(intersections$intersections@data$intersections)
    else table(intersections$intersections@data$intersections[intersections$intersections@data$intersections>0])
    x <- as.data.frame(intersections)
    x$Freq <- x$Freq / sum(x$Freq)
    names(x) <- c("intersections", "proportion")
    x$intersections <- as.integer(as.character(x$intersections))
    x$response <- study$getPrettyResponse(response)
    return(x)
  }, context=context, zeros=zeros)
  return(x)
}

intersectionsDistribution <- getIntersectionsDistribution(responses, context)
p <- ggplot(intersectionsDistribution, aes(intersections, proportion)) + geom_bar(stat="identity") +
  facet_grid(~response, scales="free_x") + xlab("Number of intersections") + ylab("Proportion") + theme_presentation(16)
print(p)
saveFigure(p, filename="IntersectionsDistribution.svg", bg="transparent")

intersectionsDistribution <- getIntersectionsDistribution(responses, context, zeros=FALSE)
p <- ggplot(intersectionsDistribution, aes(intersections, proportion)) + geom_bar(stat="identity") +
  facet_grid(~response, scales="free_x") + xlab("Number of intersections (>0)") + ylab("Proportion") + theme_presentation(16)
print(p)
saveFigure(p, filename="IntersectionsDistributionNoZero.svg", bg="transparent")

######
### CORINE habitat types
######

#colortable <- study$studyArea$habitat@legend@colortable
#colortable[colortable=="#000000"] <- "#FFFFFF"
#p <- gplot(study$studyArea$habitat, maxpixels=50000*10) + geom_raster(aes(fill=as.factor(value))) +
#  scale_fill_manual(values=colortable) +
#  coord_equal() + theme_raster()
#print(p)
#saveFigure(p, filename="CORINE.svg", dimensions=dim(study$studyArea$habitat), bg="transparent")

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
    weights$response <- study$getPrettyResponse(response)
    return(weights)
  }, context=context)
  return(x)
}

weights <- getHabitatWeights(responses=responses, context=context)
p <- ggplot(weights, aes(habitat, weights, fill=habitat)) +
  geom_bar(stat="identity") + facet_grid(~response) + scale_fill_manual(values=c("#beaed4","#ffff99","#7fc97f","#fdc086","#386cb0")) +
  xlab("") + ylab("Weight") + theme_presentation(16, axis.text.x=element_text(angle=90, hjust=1))
print(p)
saveFigure(p, filename="HabitatWeights.svg", bg="transparent")

######
### Habitat weights rasters
######

for (response in responses) {
  study <- FinlandWTCStudy$new(context=context, response=response)
  
  #habitatPreferences <- HabitatSelection$new(study=study)$loadHabitatSelection()
  #habitatWeights <- CORINEHabitatWeights$new(study=study)$setHabitatSelectionWeights(habitatPreferences)
  #rawHabitat <- raster(file.path(context$scratchDirectory, "clc2006_fi25m.tif"))
  #habitatWeights$getWeightsRaster(habitat=rawHabitat, save=TRUE)
  
  weightsRaster <- study$loadHabitatWeightsRaster()
  weightsRaster <- SpatioTemporalRaster$new(study=study, layerList=list(weightsRaster), ext="svg")
  p <- weightsRaster$plotLayer(layerName=1, plotTitle=study$getPrettyResponse(response), legendTitle="Weight")
  p <- p + theme_raster(16, legend.position=c(0.1,0.6), legend.background=element_rect(color="grey")) #, text=element_text(size=20))
  print(p)
  saveFigure(p, filename=paste("HabitatWeights-", response, ".svg", sep=""), bg="transparent")
}


######
### Distance distributions
######

distances <- list()
for (response in responses) {
  study <- FinlandWTCStudy$new(context=context, response=response)
  tracks <- study$loadTracks()
  intervals <- tracks$getSampleIntervals()
  distances <- rbind(distances, data.frame(distance=intervals$intervals$dist, logDistance=log(intervals$intervals$dist), response=study$getPrettyResponse(response)))
}

stat <- ddply(distances, .(response), summarize, mean=mean(distance), sd=sd(distance), min=min(distance), max=max(distance))
print(stat)

p <- ggplot(distances, aes(logDistance)) + geom_histogram(aes(y = ..density../sum(..density..)), binwidth=density(distances$logDistance)$bw) +
  facet_grid(~response) + xlab("log Distance (km/day)") + ylab("Proportion") + theme_presentation(16)
print(p)
saveFigure(p, filename="DistanceDistributions.svg", bg="transparent")

######
### Power law fit for simulated data
######

for (response in responses) {
  study <- FinlandWTCStudy$new(context=context, response=response)
  intervals <- study$loadSampleIntervals()
  distances <- intervals$getDistanceCurve()
  
  p <- ggplot(intervals$intervals, aes(intervalH, distanceKm)) + geom_point(alpha=.3) +
    ylab("Distance / day (km)") + xlab("Sampling interval (h)") +
    geom_line(data=distances, aes(intervalH, distanceKm), color="red", size=1) +
    geom_smooth(data=distances, aes(ymin=distanceKm25, ymax=distanceKm975), stat="identity") +
    ylim(c(0, max(intervals$intervals$distanceKm))) +
    theme_presentation() + ggtitle(study$getPrettyResponse(response))
  
  plot(p)
  saveFigure(p, filename=paste("DistanceCorrection-", response, ".svg", sep=""), bg="transparent")
  # Fixed effects: populationDensity + rrday + snow + tday
  intervals$estimatedValuesSummary()  
}

######
### Distance covariates
######

year <- 2000
covariates <- FinlandCovariates$new(study=study)
populationDensity <- covariates$loadPopulationDensityYear(year)
area <- prod(res(populationDensity) / 1000)
to <- projectExtent(populationDensity, study$studyArea$proj4string)
x <- projectRaster(populationDensity, to)

p <- gplot(x) + geom_raster(aes(fill=log(as.numeric(value))+1)) + coord_equal() +
  scale_fill_gradient(low="white", high="steelblue", na.value=NA,
                      breaks=round(seq(log(minValue(populationDensity)), log(maxValue(populationDensity)), length.out=7)),
                      guide=guide_legend(title=bquote(paste("Log density"~.(area)~km^-2)))) +
  geom_polygon(data=boundaryDF, aes(long, lat, group=group), color="black", fill=NA) +
  theme_raster(20, legend.position="right") + ggtitle(paste("Homo erectus", year))
print(p)
saveFigure(p, filename="HumanPopulationDensity.svg", bg="transparent")

weather <- covariates$loadWeatherYear(year)
weatherRaster <- SpatioTemporalRaster$new(study=study, ext="svg")
weatherRaster$addLayer(weather$month1snow, "snowcover")
weatherRaster$addLayer(weather$month1rrday, "precip")
weatherRaster$addLayer(weather$month1tday, "temp")
p <- weatherRaster$plotLayer("snowcover", boundary=boundaryDF, digits=0, plotTitle="January 2000", legendTitle="Snow depth (cm)") +
  theme_raster(20, legend.position="right")
plot(p)
saveFigure(p, filename="WeatherSnowCover.svg", bg="transparent")
p <- weatherRaster$plotLayer("precip", boundary=boundaryDF, digits=0, plotTitle="January 2000", legendTitle="Precipitation (?)") +
  theme_raster(20, legend.position="right")
plot(p)
saveFigure(p, filename="WeatherPrecipitation.svg", bg="transparent")
p <- weatherRaster$plotLayer("temp", boundary=boundaryDF, digits=0, plotTitle="January 2000", legendTitle="Temperature (°C)") +
  theme_raster(20, legend.position="right")
plot(p)
saveFigure(p, filename="WeatherTemperature.svg", bg="transparent")

######
### Predicted distance histograms and rasters
######

distances <- data.frame()
distanceRasters <- list()
for (response in responses) {
  study <- FinlandWTCStudy$new(context=context, response=response, distanceCovariatesModel=~populationDensity+rrday+snow+tday-1, trackSampleInterval=2)
  intersections <- study$loadIntersections()
  intersections$predictDistances()
  predictedDistances <- intersections$intersections$distance / 1000
    
  #intervals <- study$loadSampleIntervals()
  #estimates <- study$loadEstimates()
  #estimates$saveMeshNodeCovariates()
  #predictedDistances <- estimates$covariates$distance / 1000
  #coords <- repeatMatrix(estimates$getUnscaledMeshCoordinates(), length(estimates$years))
  
  x <- data.frame(estimates$locations, distance=predictedDistances, logDistance=log(predictedDistances), response=study$getPrettyResponse(response))
  x$year <- intersections$intersections$year
  x$year <- rep(estimates$years, each=estimates$mesh$n)
  #x <- data.frame(year=estimates$covariates$year, distance=predictedDistances, logDistance=log(predictedDistances), response=study$getPrettyResponse(response))
  distances <- rbind(distances, x)
  
  for (year in estimates$years) {
    distanceRasters[[response]] <- SpatioTemporalRaster$new(study=study, ext="svg")
    distancesYear <- x[x$year==year,]
    y <- estimates$project(distancesYear$distance)
    names(y) <- year
    distanceRasters[[response]]$rasterStack <- stack(y)
    #distanceRasters[[response]]$interpolate(distancesYear, transform=sqrt, inverseTransform=square, layerNames=year) # TODO: doesn't work well
    p <- distanceRasters[[response]]$plotLayer(paste("X", year, sep=""), boundary=boundaryDF, digits=1, plotTitle=paste(study$getPrettyResponse(response), year), legendTitle="Distance (km/day)")
    p <- p + theme_raster(20, legend.position="right")
    plot(p)
    saveFigure(p, filename=paste("DistanceRaster-", year, "-", response, ".svg", sep=""), bg="transparent")
  }
}

stat <- ddply(distances, .(response, year), summarize, mean=mean(distance), sd=sd(distance), min=min(distance), max=max(distance))
print(stat)

limits <- aes(ymax=mean+sd, ymin=mean-sd)
p <- ggplot(stat, aes(year, mean)) + geom_line(size=1) + geom_errorbar(limits, size=1) + facet_grid(~response) +
  xlab("Year") + ylab("Predicted distance (km/day), interval = 2h") + theme_presentation(16)
print(p)
saveFigure(p, filename="PredictedDistanceTimeSeries.svg")

p <- ggplot(distances, aes(logDistance)) + geom_histogram(aes(y=..density../sum(..density..)), binwidth=density(distances$logDistance)$bw) +
  facet_grid(~response, scales="free_x") + xlab("log Predicted distance (km/day), interval = 2h") + ylab("Proportion") + theme_presentation(16)
print(p)
saveFigure(p, filename="PredictedDistanceDistributions.svg")


######
### Population density
######

getPopulationDensity <- function(responses, timeModels, spatialModels, context, withHabitatWeights=FALSE) {
  populationDensity <- list()
  for (i in 1:length(responses)) {
    response <- responses[i]
    timeModel <- timeModels[i]
    spatialModel <- spatialModels[i]
    
    context <- Context$new(resultDataDirectory=wd.data.results, processedDataDirectory=wd.data.processed, rawDataDirectory=wd.data.raw, scratchDirectory=wd.scratch, figuresDirectory=wd.figures)
    study <- FinlandWTCStudy$new(context=context, response=response, distanceCovariatesModel=~populationDensity+rrday+snow+tday-1, trackSampleInterval=2)
    
    model <- if (spatialModel) FinlandSmoothModelSpatioTemporal(study=study)$setModelName("nbinomial", paste("matern", timeModel, sep="-"))
    else FinlandSmoothModelTemporal(study=study)$setModelName("nbinomial", timeModel)
    model$offsetScale <- 1000^2 # TODO: quickfix
    
    populationDensity[[response]] <- study$getPopulationDensity(model=model, withHabitatWeights=withHabitatWeights, saveDensityPlots=FALSE, getSD=FALSE)

    #study <- FinlandWTCStudy$new(context=context, response=response, distanceCovariatesModel=~populationDensity+rrday+snow+tday-1, trackSampleInterval=2)
    #populationDensity[[response]] <- study$getPopulationDensity2(withHabitatWeights=withHabitatWeights) # NOTE!!!
  }
  return(populationDensity)
}

# TODO: legend
# TODO: standard deviation, include all focal species when data available
populationDensity <- getPopulationDensity(responses=responses, timeModels=c("ar1", "ar1", "rw2"), spatialModels=c(T, T, F), context=context)
for (response in responses) {
  study <- FinlandWTCStudy$new(context=context, response=response)  
  popdens <- populationDensity[[response]]$mean
  popdens$ext="svg"
  popdens$animate(name="PopulationDensity", delay=50, ggfun=function(p) p + theme_raster(20))
}

weightedPopulationDensity <- getPopulationDensity(responses=responses, timeModels=c("ar1", "ar1", "rw2"), spatialModels=c(T, T, F), context=context, withHabitatWeights=TRUE)
for (response in responses) {
  study <- FinlandWTCStudy$new(context=context, response=response)
  wpopdens <- weightedPopulationDensity[[response]]$mean
  wpopdens$ext="svg"
  wpopdens$animate(name="WeightedPopulationDensity", delay=50, ggfun=function(p) p + theme_raster(20))
}

######
### Population size
######

populationSize <- data.frame()

for (response in responses) {
  study <- FinlandWTCStudy$new(context=context, response=response)
  names(weightedPopulationDensity[[response]]$mean$rasterStack) <- 1989:2011  # TODO: quickfix, fix the bug...
  x <- weightedPopulationDensity[[response]]$mean$integrate(volume=FinlandPopulationSize$new(study=study))
  x$loadValidationData()
  y <- x$sizeData
  y$response <- study$getPrettyResponse(response)
  print(coef(x$match()))
  if (response == "canis.lupus") y$"Adjusted estimated" <- y$Estimated * .75
  else if (response == "lynx.lynx") {
    y$"Adjusted estimated" <- y$Estimated * .3
    y$Estimated <- y$Estimated / 2
  }
  else if (response == "rangifer.tarandus.fennicus") {
    y$"Adjusted estimated" <- y$Estimated / 2
    #y$Validation <- with(y, approx(Year, Validation, n=nrow(y)))$y
  }
  populationSize <- rbind(populationSize, y)
}

x <- melt(populationSize, id.vars=c("Year","response"), variable.name="Variable")
p <- ggplot(x, aes(Year, value, group=Variable, colour=Variable, fill=Variable)) + facet_wrap(~response, scales="free_y") +
  geom_bar(subset=.(Variable=="Validation"), stat="identity") + geom_line(subset=.(Variable!="Validation"), size=2)  +
  scale_colour_manual(values=c("steelblue","violetred1","darkgreen")) +
  scale_fill_manual(values=c(NA,NA,"darkgreen")) +
  xlab("Year") + ylab("Population size") +
  theme_presentation(16, axis.text.x=element_text(angle=90, hjust=1)) + theme(legend.position="bottom")
print(p)
saveFigure(p, filename="PopulationSize.svg", bg="transparent")


######
### Straight-line distance error
######

study <- SimulationStudy$new(response="Intensive")$newInstance(context=context)
tracks <- study$loadTracks(iteration=as.integer(1))
tracks$tracks <- subset(tracks$tracks, id == 1 & year == 2001 & month == 1 & day %in% c(2,3))
intervals <- MovementSampleIntervals$new(study=study)
thinnedTracks <- intervals$getThinnedTracksSampleIntervals(tracks=tracks)

plotTracks <- data.frame()
for (i in seq(1, 11, by=5)) {
  x <- thinnedTracks$tracksList[[i]]$tracks
  x$thin <- i 
  x <- subset(x, date >= as.POSIXct("2001-01-03 00:00:00") & date <= as.POSIXct("2001-01-03 03:40:00"))
  plotTracks <- rbind(plotTracks, x)
}
plotTracks$thin <- factor(plotTracks$thin)
plotTracks$dt <- factor(plotTracks$dt / 60)

segment <- aes(x=x, y=y, xend=c(tail(x, n=-1), NA), yend=c(tail(y, n=-1), NA), group=dt, colour=dt)
ending <- arrow(length=unit(0.7, "cm"))
p <- ggplot(plotTracks, mapping=segment) +
  theme_raster(20, legend.position=c(0.5,0.85), legend.background=element_rect(color="grey")) + 
  #guides(colour=guide_legend(title=expression(paste(Delta, "t (min)")))) +
  guides(colour=guide_legend(title="Sampling interval (min)")) +
  scale_colour_manual(values=c("#4054de","#8094de","#b0c4de")) +
  ggtitle("Straight-line distance error") +
  geom_path(size=2, arrow=ending) + geom_point(size=4)
  #geom_segment(size=2, arrow=ending) # doesn't work with grouping
print(p)
saveFigure(p, filename="StraightLineDistanceError.svg", bg="transparent")


######
### Habitat usage sampling
######

usage <- HabitatSelection$new(study=study)
tracks <- study$loadTracks()
p <- usage$plotSampleSteps(tracks=tracks, index=0:3+800-69)
print(p)
saveFigure(p, filename="HabitatUsageSampling.svg", bg="transparent")


######
### Simulations
######

task_id <- 0; source(file.path(path.package("WTC"), "simulation-test", "simulate.R"))

simulate(scenario="A", nIterations=as.integer(1), plot=TRUE)
simulate(scenario="B", nIterations=as.integer(1), plot=TRUE)

simulate(scenario="D", nIterations=as.integer(1), plot=TRUE)
simulate(scenario="E", nIterations=as.integer(1), plot=TRUE)


context <- Context$new(resultDataDirectory=wd.data.results, processedDataDirectory=wd.data.processed, rawDataDirectory=wd.data.raw, scratchDirectory=wd.scratch, figuresDirectory=wd.figures)
study <- SimulationStudy$new()$newInstance(context=context, isTest=T)
png(file.path(context$figuresDirectory, "SimulationStudyArea-D-test.png"))
plot(study$studyArea$habitat)
dev.off()

p <- ggplot(study$studyArea$toGGDF(), aes(long, lat, group=group)) + geom_path(size=1) + theme_raster()
plot(p)
saveFigure(p, filename="SimulationStudyArea-A-test.svg", bg="transparent")


######
### Russian data
######

# TODO: awful code, fix!

response <- "canis.lupus"
study <- RussiaWTCStudy(context=context, response=response)
study$studyArea$saveBoundary(study=study)
load(file=study$studyArea$getBoundaryFileName(), envir=as.environment(study$studyArea))
plot(study$studyArea$boundary)

intersections <- study$loadIntersections()
minx <- xmin(extent(study$studyArea$habitat))
maxx <- xmax(extent(study$studyArea$habitat))
miny <- ymin(extent(study$studyArea$habitat))
maxy <- ymax(extent(study$studyArea$habitat))
window <- matrix(c(minx,miny, minx,maxy, maxx,maxy, maxx,miny, minx,miny), ncol=2, byrow=T)
coords <- coordinates(intersections$intersections)
z <- intersections$intersections[point.in.polygon(coords[,1], coords[,2], window[,1], window[,2])==1,]
plot(z)
boundary <- subset(study$studyArea$boundary, NAME_LAT %in% z$District_Lat & ADM4_ID %in% z$RegionID)
plot(boundary)

x <- subset(z, Year==2001)
x$density <- x$intersections / (x$length/1000) / x$duration #/ x$area # / track length
x$NAME_LAT <- x$District_Lat
x$ADM4_ID <- x$RegionID
y <- boundary
ids <- unlist(lapply(y@polygons, function(x) { x@ID }))
y$id <- ids
# TODO: merge so that missing districts have missing values
y@data <- merge(y@data, x@data[,c("NAME_LAT","ADM4_ID","density")], by=c("NAME_LAT", "ADM4_ID"))
y@polygons <- y@polygons[ids %in% y$id]
#plot(y, col=colorRampPalette(c('white', 'red'))(length(x$density))[rank(x$density)])
#legend("right", "(x,y)", title="Wolf crossing density")


library(ggplot2)
library(scales)
library(plyr)
#spChFIDs(y, as.character(1:length(y)))
#y$id <- unlist(lapply(y@polygons, function(x) { x@ID }))

#y$density[y$density>8e-5] <- NA
yf <- fortify(y)
yfd <- join(yf, y@data[,c("id","density")], by="id")
p <- ggplot(yfd, aes(long, lat, group=group)) +
  geom_polygon(aes(fill=density)) + scale_fill_continuous("Crossing density / km", low="white", high=muted("red")) +
  geom_path(color="gray", linestyle=2) +
  theme_raster(20, legend.position="bottom", plot.margin=unit(c(0,0,-1,2), "lines"),
               legend.key=element_rect(size=0.5), 
               legend.key.height=unit(1, "cm"), 
               legend.key.width=unit(1, "cm")) +
  ggtitle("Волк")
saveFigure(p, filename="CrossingDensity-canis.lupus-Russia.svg", bg="transparent")
print(p)
