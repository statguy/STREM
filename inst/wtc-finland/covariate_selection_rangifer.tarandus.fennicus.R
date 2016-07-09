library(parallel)
library(doMC)
registerDoMC(cores=detectCores())
library(STREM)
source("~/git/STREM/setup/WTC-Boot.R")
library(SpaceTimeModels)
library(stringr)

#response <- "canis.lupus"
#response <- "lynx.lynx"
response <- "rangifer.tarandus.fennicus"


context <- Context(resultDataDirectory=wd.data.results, processedDataDirectory=wd.data.processed, rawDataDirectory=wd.data.raw, scratchDirectory=wd.scratch, figuresDirectory=wd.figures)
study <- FinlandWTCStudy(context=context, response=response)
intersections <- FinlandWTCIntersections$new(study=study)$loadIntersections()

# Quickfix for the errorneously saved covariate data
study2 <- FinlandWTCStudy(context=context, response="canis.lupus")
x <- FinlandWTCIntersections$new(study=study2)$setCovariatesId(tag="density")$loadCovariates()
y <- plyr::join(as.data.frame(intersections$intersections), as.data.frame(x$intersections), by=c("x","y","id","year"))
intersections$intersections <- SpatialPointsDataFrame(y[,c("x","y")], data=y[,!names(y) %in% c("y","x")], proj4string=CRS("+init=epsg:2393"))
names(intersections$intersections) <- stringr::str_replace_all(names(intersections$intersections), "[\\[\\]_]", ".")
intersections$covariateNames <- stringr::str_subset(names(intersections$intersections), "^(habitat|topography|human)")
intersections$removeMissingCovariates()

# Use durations of 1 only
index <- intersections$intersections$duration == 1
intersections$intersections <- intersections$intersections[index,]
nrow(intersections$intersections@data)

# Compress habitat covariates
dataNames <- colnames(intersections$intersections@data)
for (scale in 2^(0:10) * 62.5) {
  variables <- stringr::str_subset(dataNames, paste0("^habitat\\.\\w+\\.\\.", scale, "$"))
  intersections$compressCovariates(variables, name=paste0("habitat.", scale))
}

# Scale covariates
variables <- c(stringr::str_subset(dataNames, "^topography\\.rel\\.\\d+$"), stringr::str_subset(dataNames, "^human\\.\\d+$"))
intersections$scaleCovariates(variables)
names(intersections$intersections@data)

# Convert spatio-temporal data to spatial
data <- ddply(cbind(coordinates(intersections$intersections), intersections$intersections@data), .(id),
              function(x) data.frame(x=x[1,1], y=x[1,2], length=sum(x$length), distance=sum(x$distance), intersections=sum(x$intersections), x[1,12:ncol(x)]))

mesh <- NonConvexHullMesh$new(SpatialPoints(data[,c("x","y")]), knotsScale=1e6, cutoff=0.01e6, offset=c(1,0.1)*1e6, maxEdge=c(0.05,1)*1e6, convex=0.12)
mesh$getINLAMesh()$n
tracks <- study$loadTracks()
distance <- mean(tracks$getDistances())


results <- list()

smoothOnly <- function(response, mesh, data, distance, formula) {
  model <- ContinuousSpaceModel$new(offsetScale=1e7)
  model$setSpatialMesh(mesh)
  model$setSpatialPrior()
  
  model$setSmoothingModel()
  model$addObservationStack(sp=SpatialPoints(data[,c("x","y")]), response=data$intersections, offset=data$length * distance, covariates=data, tag="obs")
  model$setLikelihood("nbinomial")
  model$estimate(verbose=F)
  print(summary(model$getResult()))
  
  x <- data.frame(response=response, formula="", waic=model$getResult()$waic$waic)
  results <<- rbind(results, x)
  return(invisible(model))
}

smoothOnly(response, mesh, data, distance, "")
# 570.15

selectCovariates <- function(response, mesh, data, distance, formula) {
  #distance <- data$distance
  #distance <- mean(intersections$intersections@data$distance)
  #sdDistance <- sd(intersections$intersections@data$distance)
  model <- ContinuousSpaceModel$new(offsetScale=1e7)
  model$setSpatialMesh(mesh)
  model$setSpatialPrior()
  
  #model$setSmoothingModel()
  #model$addObservationStack(sp=SpatialPoints(data[,c("x","y")]), response=data$intersections, offset=data$length * distance, covariates=data, tag="obs")
  model$setLikelihood("nbinomial")
  #model$estimate(verbose=F)
  #summary(model$getResult())
  
  names(data)
  model$setCovariatesModel(formula, data)
  model$addObservationStack(sp=SpatialPoints(data[,c("x","y")]), response=data$intersections, offset=data$length * distance, covariates=data, tag="obs")
  model$estimate(verbose=F)
  print(summary(model$getResult()))
  
  x <- data.frame(response=response, formula=as.character(formula)[2], waic=model$getResult()$waic$waic)
  results <<- rbind(results, x)
  return(invisible(model))
}

selectCovariates(response, mesh, data, distance, ~PC1.habitat.62.5+PC2.habitat.62.5+PC1.habitat.125+PC2.habitat.125+PC1.habitat.250+PC2.habitat.250+PC1.habitat.500+PC2.habitat.500+PC1.habitat.1000+PC2.habitat.1000+PC1.habitat.2000+PC2.habitat.2000+PC1.habitat.4000+PC2.habitat.4000+PC1.habitat.8000+PC2.habitat.8000+PC1.habitat.16000+PC1.habitat.32000+PC1.habitat.64000+topography.rel.1000+topography.rel.2000+topography.rel.4000+topography.rel.8000+topography.rel.16000+topography.rel.32000+topography.rel.64000+human.1000+human.2000+human.4000+human.8000+human.16000+human.32000+human.64000)
# 669.15

result <- selectCovariates(response, mesh, data, distance, ~habitat.Agriculture..8000+habitat.Urban..250+habitat.Forestland..32000+topography.rel.32000)
# 658.14
result <- selectCovariates(response, mesh, data, distance, ~habitat.Urban..250+habitat.Forestland..32000+topography.rel.32000)
