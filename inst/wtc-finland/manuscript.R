#+ setup, include=FALSE
library(knitr)
opts_chunk$set(echo=FALSE, warning=FALSE, message=FALSE, tidy=FALSE, fig.width=10, fig.height=15, cache=TRUE)

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
###
######

