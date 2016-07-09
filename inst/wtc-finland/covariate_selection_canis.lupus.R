library(parallel)
library(doMC)
registerDoMC(cores=detectCores())
library(STREM)
source("~/git/STREM/setup/WTC-Boot.R")
library(SpaceTimeModels)
library(stringr)

response <- "canis.lupus"

context <- Context(resultDataDirectory=wd.data.results, processedDataDirectory=wd.data.processed, rawDataDirectory=wd.data.raw, scratchDirectory=wd.scratch, figuresDirectory=wd.figures)
study <- FinlandWTCStudy(context=context, response=response)
intersections <- FinlandWTCIntersections$new(study=study)$loadIntersections()

# Quickfix for the errorneously saved covariate data
study2 <- FinlandWTCStudy(context=context, response="canis.lupus")
x <- FinlandWTCIntersections$new(study=study2)$setCovariatesId(tag="density")$loadCovariates()
y <- plyr::join(as.data.frame(intersections$intersections), as.data.frame(x$intersections))
intersections$intersections <- SpatialPointsDataFrame(y[,c("x","y")], data=y[,!names(y) %in% c("y","x")], proj4string=CRS("+init=epsg:2393"))
names(intersections$intersections) <- stringr::str_replace_all(names(intersections$intersections), "[\\[\\]_]", ".")
intersections$covariateNames <- stringr::str_subset(names(intersections$intersections), "^(habitat|topography|human)")
intersections$removeMissingCovariates()

index <- coordinates(intersections$intersections)[,2] < 7.3e6 # Cut off Lapland since all wolves are killed there
index <- index & intersections$intersections$duration == 1 # Use durations of 1 only
intersections$intersections <- intersections$intersections[index,]

# Scale covariates
variables <- c(stringr::str_subset(dataNames, "^topography\\.rel\\.\\d+$"), stringr::str_subset(dataNames, "^human\\.\\d+$"))
intersections$scaleCovariates(variables)
names(intersections$intersections@data)

# 0.3 + 0.7 = 1
# 0.2 + 0.8 = 1
# 0.3 / (0.3+0.2) = 0.6, 0.7/(0.7+0.8) = 0.4666 -- 0.6 + 0.466 !=
# 0.2 / (0.3+0.2) = 0.4, 0.8/(0.7+0.8) = 0.5333 -- 0.4 + 0.533 !=

# J = variables
# S = scales

# J * S = scaled variables
# J = new variables from PCA


# Compress covariates
dataNames <- colnames(intersections$intersections@data)
for (scale in 2^(0:3) * 62.5) {
  variables <- stringr::str_subset(dataNames, paste0("^habitat\\.\\w+\\.\\.", scale, "$"))
  intersections$compressCovariates(variables, prefix="PC", name=paste0("habitat.", scale))
}
for (scale in 2^(4:10) * 62.5) {
  variables <- stringr::str_subset(dataNames, paste0("^(habitat\\.\\w+\\.\\.", scale, "|topography\\.rel\\.", scale, "\\.scaled|human\\.", scale, "\\.scaled)$"))
  intersections$compressCovariates(variables, prefix="PC", name=paste0("all.", scale))
}
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
# 1074.72

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

# Method 1: use PC variables
selectCovariates(response, mesh, data, distance, reformulate(stringr::str_subset(colnames(data), "^PC\\.")))
# 1074.78
selectCovariates(response, mesh, data, distance, ~PC.PC1.habitat.250+PC.PC2.all.64000)
# 1072.57
selectCovariates(response, mesh, data, distance, ~PC.PC1.habitat.250)
# 1070.76


# Method 2: fix baseline scale
selectCovariates(response, mesh, data, distance, ~habitat.Urban..62.5+habitat.Agriculture..62.5+habitat.Forestland..62.5+habitat.Peatland..62.5+topography.rel.1000+human.1000)
# 1079.32
selectCovariates(response, mesh, data, distance, ~habitat.Urban..125+habitat.Agriculture..125+habitat.Forestland..125+habitat.Peatland..125+topography.rel.1000+human.1000)
# 1078.34
selectCovariates(response, mesh, data, distance, ~habitat.Urban..250+habitat.Agriculture..250+habitat.Forestland..250+habitat.Peatland..250+topography.rel.1000+human.1000)
# 1077.23
selectCovariates(response, mesh, data, distance, ~habitat.Urban..500+habitat.Agriculture..500+habitat.Forestland..500+habitat.Peatland..500+topography.rel.1000+human.1000)
# 1078.71
selectCovariates(response, mesh, data, distance, ~habitat.Urban..1000+habitat.Agriculture..1000+habitat.Forestland..1000+habitat.Peatland..1000+topography.rel.1000+human.1000)
# 1081.90
selectCovariates(response, mesh, data, distance, ~habitat.Urban..2000+habitat.Agriculture..2000+habitat.Forestland..2000+habitat.Peatland..2000+topography.rel.2000+human.2000)
# 1083.05
selectCovariates(response, mesh, data, distance, ~habitat.Urban..4000+habitat.Agriculture..4000+habitat.Forestland..4000+habitat.Peatland..4000+topography.rel.4000+human.4000)
# 1079.92
selectCovariates(response, mesh, data, distance, ~habitat.Urban..8000+habitat.Agriculture..8000+habitat.Forestland..8000+habitat.Peatland..8000+topography.rel.8000+human.8000)
# 1076.78
selectCovariates(response, mesh, data, distance, ~habitat.Urban..16000+habitat.Agriculture..16000+habitat.Forestland..16000+habitat.Peatland..16000+topography.rel.16000+human.16000)
# 1067.00
selectCovariates(response, mesh, data, distance, ~habitat.Urban..32000+habitat.Agriculture..32000+habitat.Forestland..32000+habitat.Peatland..32000+topography.rel.32000+human.32000)
# 1061.45
selectCovariates(response, mesh, data, distance, ~habitat.Urban..64000+habitat.Agriculture..64000+habitat.Forestland..64000+habitat.Peatland..64000+topography.rel.64000+human.64000)
# 1073.94



# Method 3: guess
result <- selectCovariates(response, mesh, data, distance, ~habitat.Agriculture..8000+habitat.Urban..250+habitat.Forestland..32000+topography.rel.32000)
# 1060.09

selectCovariates(response, mesh, data, distance, ~PC1.habitat.62.5+PC2.habitat.62.5+PC1.habitat.125+PC2.habitat.125+PC1.habitat.250+PC2.habitat.250+PC1.habitat.500+PC2.habitat.500+PC1.habitat.1000+PC2.habitat.1000+PC1.habitat.2000+PC2.habitat.2000+PC1.habitat.4000+PC2.habitat.4000+PC1.habitat.8000+PC2.habitat.8000+PC1.habitat.16000+PC1.habitat.32000+PC1.habitat.64000+topography.rel.1000+topography.rel.2000+topography.rel.4000+topography.rel.8000+topography.rel.16000+topography.rel.32000+topography.rel.64000+human.1000+human.2000+human.4000+human.8000+human.16000+human.32000+human.64000)
#
selectCovariates(response, mesh, data, distance, ~PC1.habitat.250+PC2.habitat.250+topography.rel.2000+topography.rel.4000+topography.rel.8000+topography.rel.32000+human.1000)
#
selectCovariates(response, mesh, data, distance, ~topography.rel.2000+topography.rel.4000+topography.rel.8000+topography.rel.32000)
#
selectCovariates(response, mesh, data, distance, ~topography.rel.2000+topography.rel.4000+topography.rel.8000+topography.rel.32000 +habitat.Urban..250)
#
selectCovariates(response, mesh, data, distance, ~topography.rel.2000+topography.rel.4000+topography.rel.8000+topography.rel.32000 +habitat.Urban..250+habitat.Agriculture..8000)
#
selectCovariates(response, mesh, data, distance, ~topography.rel.2000+topography.rel.4000+topography.rel.8000+topography.rel.32000 +habitat.Urban..250+habitat.Agriculture..8000+habitat.Forestland..32000)
#
selectCovariates(response, mesh, data, distance, ~topography.rel.2000+topography.rel.4000+topography.rel.8000+topography.rel.32000 +habitat.Urban..250+habitat.Agriculture..8000+habitat.Forestland..32000+habitat.Peatland..250)
# 1054.47
selectCovariates(response, mesh, data, distance, ~topography.rel.4000+topography.rel.8000+topography.rel.32000 +habitat.Urban..250+habitat.Agriculture..8000+habitat.Forestland..32000+habitat.Peatland..250)
# 1059.59
selectCovariates(response, mesh, data, distance, ~topography.rel.2000+topography.rel.32000 +habitat.Urban..250+habitat.Agriculture..8000+habitat.Forestland..32000+habitat.Peatland..250)
# 1061.40
selectCovariates(response, mesh, data, distance, ~topography.rel.2000+topography.rel.4000+topography.rel.8000+topography.rel.32000 +habitat.Urban..250+habitat.Agriculture..8000+habitat.Forestland..32000+habitat.Peatland..250 + human.1000)
# 1057.94
selectCovariates(response, mesh, data, distance, ~topography.rel.2000+topography.rel.4000+topography.rel.8000+topography.rel.32000 +habitat.Urban..250+habitat.Agriculture..8000+habitat.Forestland..32000+habitat.Peatland..250 + human.32000)
# 1054.12
selectCovariates(response, mesh, data, distance, ~topography.rel.2000+topography.rel.4000+topography.rel.8000+topography.rel.32000 +habitat.Urban..250+habitat.Agriculture..8000+habitat.Forestland..2000+habitat.Peatland..250 + human.32000)
# 1050.34
selectCovariates(response, mesh, data, distance, ~topography.rel.2000+topography.rel.4000+topography.rel.8000+topography.rel.32000 +habitat.Urban..250+habitat.Agriculture..1000+habitat.Forestland..2000+habitat.Peatland..250 + human.32000)
# 1054.36
selectCovariates(response, mesh, data, distance, ~topography.rel.2000+topography.rel.4000+topography.rel.8000+topography.rel.32000 +habitat.Urban..125+habitat.Agriculture..8000+habitat.Forestland..2000+habitat.Peatland..250 + human.32000)
# 1050.32
selectCovariates(response, mesh, data, distance, ~topography.rel.2000+topography.rel.4000+topography.rel.8000+topography.rel.32000 +habitat.Urban..125+habitat.Agriculture..32000+habitat.Forestland..2000+habitat.Peatland..250 + human.32000)
# 1053.50
selectCovariates(response, mesh, data, distance, ~topography.rel.1000+topography.rel.4000+topography.rel.8000+topography.rel.32000 +habitat.Urban..125+habitat.Agriculture..32000+habitat.Forestland..2000+habitat.Peatland..250 + human.32000)
# 1057.99
selectCovariates(response, mesh, data, distance, ~topography.rel.1000+topography.rel.4000+topography.rel.8000+topography.rel.32000 +habitat.Urban..125+habitat.Agriculture..32000+habitat.Forestland..2000+habitat.Forestland..32000+habitat.Peatland..250 + human.32000)
# 1061.39






cor(data$topography.rel.2000, data$habitat.Peatland..1000)
cor(data$habitat.Urban..1000, data$human.1000)

cbind(data$intersections, round(result$getFittedResponse()$responseMean), round(result2$getFittedResponse()$responseMean))
cor(data$intersections, result$getFittedResponse()$responseMean)
cor(data$intersections, result2$getFittedResponse()$responseMean)

result$getFittedLinearPredictor()$etaMean
result2$getFittedLinearPredictor()$etaMean
