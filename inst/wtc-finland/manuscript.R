#+ setup, include=FALSE
library(knitr)
opts_chunk$set(echo=FALSE, warning=FALSE, message=FALSE, tidy=FALSE, fig.width=10, fig.height=15, cache=TRUE)#, cache.path="~/tmp/cache")

#' NOTE: Figures and tables may be out of sync with the manuscript!
#' ----------------------------------------------------------------

```{r initialize}
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
study$studyArea$loadBoundary(thin=T, tolerance=0.001)
boundaryDF <- ggplot2::fortify(study$studyArea$boundary)
responses <- c("canis.lupus", "lynx.lynx", "rangifer.tarandus.fennicus")
```

######
### Survey routes (Figure 1)
######

```{r survey-routes}
surveyRoutes <- study$loadSurveyRoutes(findLengths=FALSE)
surveyRoutesSP <- SpatialLinesDataFrame(surveyRoutes$surveyRoutes, data=data.frame(x=1:length(surveyRoutes$surveyRoutes)), match.ID=FALSE)
surveyRoutesDF <- ggplot2::fortify(surveyRoutesSP)
p <- ggplot(boundaryDF, aes(long, lat, group=group)) + geom_polygon(colour="darkgray", fill=NA) +
  geom_polygon(data=surveyRoutesDF, aes(long, lat, group=group), colour="black", fill=NA) +
  coord_equal() + theme_raster()
print(p)
#saveFigure(p, filename="SurveyRoutes.svg", bg="transparent")
```

######
### Number of counted survey routes
######

```{r number-of-counted-survey-routes}
intersections <- study$loadIntersections(predictDistances=FALSE)
x <- ddply(as.data.frame(intersections$intersections), .(year), nrow)
s <- sd(x$V1)
m <- mean(x$V1)
cat("mean = ", m, "\n")
cat("sd = ", s, "\n")
cat("range = ", m+c(-s,s), "\n")
```

######
### % of intersections D > 1
######

```{r percent-of-intersections-day-more-than-one, message=TRUE}
intersections <- FinlandWTCIntersections$new(study=study, maxDuration=1)
intersections$saveIntersections()
```

######
### Crossing rates (D > 1 removed) (Figure 2)
######

```{r crossing-rates, fig.height=5}
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
  xlab("Year") + ylab("Crossings rate") + theme_presentation(20)
print(p)
#saveFigure(p, filename="IntersectionsRate.svg", bg="transparent")
```

######
### GPS figures
######

```{r gps-in-numbers, message=TRUE}
getTracks <- function(responses, context) {
  library(plyr)
  library(ggplot2)
  l_ply(responses, function(response, context) {
    study <- FinlandWTCStudy$new(context=context, response=response)
    tracks <- FinlandWTCTracks$new(study=study)
    tracks$loadTracks()
    
    message("Response = ", response)
    message("Number of individuals in the cleaned data = ", length(unique(id(tracks$tracks))))
    d <- as.POSIXlt(ld(tracks$tracks)$date)
    message("Years in the cleaned data = ", paste(sort(unique(d$year+1900)), collapse=","), " (", length(unique(d$year)), ")")
    message("Number of samples in the cleaned data = ", length(d))
    message("Number of days in the cleaned data = ", length(unique(paste(d$year, d$yday))))    
    
  }, context=context)
}

tracks <- getTracks(responses=responses, context=context)
```

######
### GPS movements (Figure 3)
######

```{r gps-movements, fig.height=5}
getTracks <- function(responses, context) {
  library(plyr)
  library(ggplot2)
  x <- ldply(responses, function(response, context) {
    study <- FinlandWTCStudy$new(context=context, response=response)
    tracks <- FinlandWTCTracks$new(study=study)
    tracks$loadTracks()
    tracksSP <- tracks$getSpatialLines(variables=.(id,year,burst))
    
    x <- tracksSP@data
    if (!"sex" %in% colnames(x)) x$sex <- NA
    x$juvenile <- if ("born" %in% colnames(x)) factor(ifelse(x$year - x$born <= 1, "juvenile", "adult")) else factor(NA)
    x$id <- sapply(tracksSP@lines, function(x) x@ID)
    x$sexjuv <- factor(paste(x$sex, x$juvenile, sep=", "))
    tracksDF <- ggplot2::fortify(tracksSP)
    tracksDF <- plyr::join(tracksDF, x, by="id")
    tracksDF$response <- study$getPrettyResponse(response)
    
    return(tracksDF)
  }, context=context)
  return(x)
}

tracks <- getTracks(responses=responses, context=context)
#levels(tracks$sexjuv)

p <- ggplot(tracks, aes(long, lat, group=group)) +
  geom_polygon(data=boundaryDF, aes(long, lat, group=group), colour="darkgray", fill=NA) +
  geom_path(aes(colour=sexjuv, linetype=sexjuv), alpha=.5) +
  scale_color_manual("Sex, age", values=c("darkred","red","red","darkblue","blue","blue","black")) + scale_linetype_manual("Sex, age", values=c("solid","dashed","solid","solid","dashed","solid","solid")) +
  facet_grid(~response) + theme_raster(18, legend.position="bottom", aspect.ratio=1) + coord_cartesian(ylim=range(tracks$lat)+c(-1,1)*50000)
print(p)
#saveFigure(p, filename="GPS.svg", bg="transparent")
```

######
### Straight-line-distance-error (Figure 4)
######

```{r straight-line-distance-error, fig.height=5}
ss <- SimulationStudy$new(response="Intensive")$newInstance(context=context)
tracks <- ss$loadTracks(iteration=as.integer(1))
tracks$tracks <- subset(tracks$tracks, id == 1 & year == 2001 & month == 1 & day %in% c(2,3))
intervals <- MovementSampleIntervals$new(study=ss)
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
  scale_colour_manual(values=c("gray70","gray40","gray10")) +
  #ggtitle("Straight-line distance error") +
  geom_path(size=2, arrow=ending) + geom_point(size=4)
#geom_segment(size=2, arrow=ending) # doesn't work with grouping
print(p)
#saveFigure(p, filename="StraightLineDistanceError.svg", bg="transparent")
```

######
### Raw GPS sample intervals
######

```{r raw-gps-sample-intervals}
intervals <- llply(responses, function(response, context) {
  study <- FinlandWTCStudy$new(context=context, response=response)
  tracks <- FinlandWTCTracks$new(study=study)
  tracks$loadTracks()
  intervals <- MovementSampleIntervals$new(study=study)
  intervals$loadSampleIntervals()
  return(table(intervals$intervals$intervalH))
}, context=context)
intervals
```

######
### Power law fits for distance - sampling interval (Figure 5)
######

```{r distance-fits, fig.height=5}
distances <- ldply(responses, function(response, context) {
  study <- FinlandWTCStudy$new(context=context, response=response)
  intervals <- study$loadSampleIntervals()
  x <- intervals$intervals
  x$response <- study$getPrettyResponse(response=response)
  return(x)
}, context=context)

intervals <- ldply(responses, function(response, context) {
  study <- FinlandWTCStudy$new(context=context, response=response)
  intervals <- study$loadSampleIntervals()
  distances <- intervals$getDistanceCurve()
  distances$response <- study$getPrettyResponse(response=response)
  return(distances)
}, context=context)

p <- ggplot(distances, aes(intervalH, distanceKm)) + geom_point(alpha=.3) +
  ylab("Distance / day (km)") + xlab("Sampling interval (h)") +
  geom_line(data=intervals, aes(intervalH, distanceKm), color="black", size=1) +
  geom_smooth(data=intervals, aes(ymin=distanceKm25, ymax=distanceKm975), stat="identity") +
  ylim(c(0, max(distances$distanceKm))) +
  facet_grid(~response) + theme_presentation()
print(p)
#saveFigure(p, filename=paste("DistanceCorrection-", response, ".svg", sep=""), bg="transparent")
```

######
### Observed distances - estimated distances (Table 1)
######

```{r distance-estimates-init, results="hide"}
observedDistances <- ldply(responses, function(response, context) {
  study <- FinlandWTCStudy$new(context=context, response=response)
  tracks <- study$loadTracks()
  intervals <- tracks$getSampleIntervals()
  return(data.frame(distance=intervals$intervals$dist, logDistance=log(intervals$intervals$dist), response=study$getPrettyResponse(response)))
}, context=context)

predictedDistances <- ldply(responses, function(response, context) {
  study <- FinlandWTCStudy$new(context=context, response=response, distanceCovariatesModel=~populationDensity+rrday+snow+tday-1, trackSampleInterval=2)
  intersections <- study$loadIntersections()
  intersections$predictDistances()
  distances <- intersections$intersections$distance / 1000
  return(data.frame(distance=distances, logDistance=log(distances), response=study$getPrettyResponse(response)))  
}, context=context)
```

```{r distance-estimates}
cat("Observed distances (km/day):\n")
ddply(observedDistances, .(response), summarize, mean=mean(distance), sd=sd(distance), min=min(distance), max=max(distance))
cat("Predicted distances (km/day):\n")
ddply(predictedDistances, .(response), summarize, mean=mean(distance), sd=sd(distance), min=min(distance), max=max(distance))
cat("Prediction interval = 2h")
```

# Distances from literature:
# Wolf: 7-25 km/day
# Lynx:
# FFR: 


######
### Power law fit estimates for distance - sampling interval (Table 1)
######
# Fixed effects: human population density, precipitation, snow cover, temperature

```{r distance-fit-estimates}
# Fixed effects: populationDensity + rrday + snow + tday
intervals <- l_ply(responses, function(response, context) {
  study <- FinlandWTCStudy$new(context=context, response=response)
  intervals <- study$loadSampleIntervals()
  cat("response = ", study$getPrettyResponse(response=response), "\n")
  intervals$estimatedValuesSummary()
  cat("\n\n")
}, context=context)
```

######
### Crude distance estimates (Table 1)
#####

```{r distance-crude-estimates}
responses <- c("canis.lupus","lynx.lynx","rangifer.tarandus.fennicus")
for (response in responses) {
  context <- Context$new(resultDataDirectory=wd.data.results, processedDataDirectory=wd.data.processed, rawDataDirectory=wd.data.raw, scratchDirectory=wd.scratch, figuresDirectory=wd.figures)
  study <- FinlandWTCStudy$new(context=context, response=response, distanceCovariatesModel=~populationDensity+rrday+snow+tday-1, trackSampleInterval=2)
  
  populationSize <- FinlandPopulationSize$new(study=study)
  populationSize$loadValidationData()
  
  intersections <- study$loadIntersections()
  
  A <- 338424 # km^2
  M <- ddply(intersections$getData(), .(year), function(x) sum(x$length))$V1 / 1000
  N <- populationSize$sizeData$Validation
  x <- ddply(intersections$getData(), .(year), function(x) sum(x$intersections))$V1
  
  C <- 2/pi * A / M
  l <- C * x / N
  cat(response, ": crude distance mean = ", mean(l, na.rm=T), " sd = ", sd(l, na.rm=T), " km/d\n")
    
  h <- study$loadHabitatWeightsRaster()
  h.masked <- mask(h, study$studyArea$boundary)
  n.cell <- length(h.masked) - freq(h.masked, value=NA)
  A.discrete <- n.cell * prod(res(h.masked)) / 1000^2
  h.cell <- cellStats(h.masked, sum)
  A.h <- h.cell * prod(res(h.masked)) / 1000^2 * A / A.discrete # A.discrete is a little bit smaller than A due to discretization, so apply a correction
  
  C <- 2/pi * A.h / M
  l <- C * x / N
  cat(response, ": crude weighted distance mean = ", mean(l, na.rm=T), " sd = ", sd(l, na.rm=T), " km/d\n")
  
  cat(response, ": observed number of intersections:\n")
  print(1/C * l)
}
```

######
### Habitat weights (Figure 6)
######

```{r habitat-weights, fig.height=5, results="hide"}
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
  geom_bar(stat="identity") + facet_grid(~response) + scale_fill_manual(values=c("gray40","gray40","gray40","gray40","gray40")) +
  xlab("") + ylab("Weight") + theme_presentation(20, axis.text.x=element_text(angle=90, hjust=1))
print(p)
#saveFigure(p, filename="HabitatWeights.svg", bg="transparent")
```

######
### Total population size estimates (Figure 7)
######

```{r population-size-init, results="hide"}
getPopulationDensity <- function(responses, context, withHabitatWeights=FALSE) {
  populationDensity <- list()
  for (response in responses) {
    #study <- FinlandWTCStudy$new(context=context, response=responses[1], distanceCovariatesModel=~populationDensity+rrday+snow+tday-1, trackSampleInterval=2)
    #estimates <- study$loadEstimates()
    #estimates$offsetScale <- 1000^2
    #estimates$collectEstimates()
    #estimates$saveEstimates()
    
    study <- FinlandWTCStudy$new(context=context, response=response, distanceCovariatesModel=~populationDensity+rrday+snow+tday-1, trackSampleInterval=2)
    populationDensity[[response]] <- study$getPopulationDensity(withHabitatWeights=withHabitatWeights)
  }
  return(populationDensity)
}

weightedPopulationDensity <- getPopulationDensity(responses=responses, context=context, withHabitatWeights=TRUE)
```

```{r population-size, fig.height=5}
populationSize <- data.frame()
for (response in responses) {
  study <- FinlandWTCStudy$new(context=context, response=response)
  x <- weightedPopulationDensity[[response]]$mean$integrate(volume=FinlandPopulationSize$new(study=study))
  x$loadValidationData()
  y <- x$sizeData
  y$response <- study$getPrettyResponse(response)
  print(coef(x$match()))
  if (response == "canis.lupus") y$"Transformed estimated" <- y$Estimated * .3
  else if (response == "lynx.lynx") y$"Transformed estimated" <- y$Estimated * .2
  else if (response == "rangifer.tarandus.fennicus") y$"Transformed estimated" <- y$Estimated * 1.27514e-132
  
  if (response == "rangifer.tarandus.fennicus") {
    y$Estimated <- NA
    y$Validation <- with(y, approx(Year, Validation, n=nrow(y)))$y
  }
  
  populationSize <- rbind(populationSize, y)
}

x <- melt(populationSize, id.vars=c("Year","response"), variable.name="Variable")
p <- ggplot(x, aes(Year, value, group=Variable, colour=Variable)) + geom_line(size=1) + facet_wrap(~response, scales="free_y") +
  scale_colour_manual(values=c("steelblue","violetred1","steelblue1")) +
  xlab("Year") + ylab("Population size") +
  theme_presentation(20, axis.text.x=element_text(angle=90, hjust=1)) + theme(legend.position="bottom")
print(p)
#saveFigure(p, filename="PopulationSize.svg", bg="transparent")
```
