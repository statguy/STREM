repeatMatrix <- function(m, times) do.call("rbind", replicate(times, m, simplify=FALSE))

# Converts ltraj object to SpatialPointsDataFrame object without stripping information.
# Allows setting projection as well.
ltraj2spdf2 <- function(ltraj, proj4string="") {
  library(adehabitatLT)
  library(sp)
  df <- ld(ltraj)
  coordinates(df) <- ~x+y
  proj4string(df) <- proj4string
  return(df)
}

breakDownDate <- function(x) {
  x <- as.POSIXlt(x)
  return(data.frame(year=x$year+1900, month=x$mon+1, day=x$mday, yday=x$yday))
}

euclidean <- function(c1, c2) sqrt((c1[,1] - c2[,1])^2 + (c1[,2] - c2[,2])^2)

square <- function(x) x^2

rotate <- function(x, y, angle) {
  return(cbind(x * cos(angle) - y * sin(angle), x * sin(angle) + y * cos(angle)))
}

getVector <- function(coords, distance, angle) {
  return(cbind(coords[,1] + distance * cos(angle), coords[,2] + distance * sin(angle)))
}

getTriangles <- function(centroids, angles, sideLength) {
  library(plyr)
  library(sp)
  
  l2 <- sideLength / 2
  lmid <- l2 * sin(pi / 3)
  xy <- coordinates(centroids)
  c1 <- rotate(-l2, lmid, angles) + xy
  c2 <- rotate(l2, lmid, angles) + xy
  c3 <- rotate(0, -lmid, angles) + xy
  
  triangles <- llply(1:nrow(c1), function(i, c1, c2, c3) {
    return(Lines(list(Polygon(rbind(c1[i,], c2[i,], c3[i,], c1[i,]))), ID=i))
  }, c1=c1, c2=c2, c3=c3)
  triangles <- SpatialLines(triangles, proj4string=centroids@proj4string)
  
  return(triangles)
}

getINLAEstimates <- function(marginal, fun=identity, coordsScale=1) {
  m <- inla.tmarginal(function(x) fun(x) / coordsScale, marginal)
  e <- inla.emarginal(identity, m)
  e2 <- inla.emarginal(function(x) x^2, m)
  sd <- sqrt(e2 - e^2)
  q <- inla.qmarginal(c(0.025, 0.5, 0.975), m)
  mode <- inla.mmarginal(m)
  x <- data.frame(e=e, sd=sd, q1=q[1], q2=q[2], q3=q[3], mode=mode)
  colnames(x) <- c("mean", "sd", "0.025quant","0.5quant","0.975quant", "mode")
  return(x)
}

inverseDistanceWeightningImpute <- function(data, varname, formula=as.formula(paste(varname, "~1")), locations=~x+y) {
  library(gstat)
  data.full <- data; data.full <- data.full[complete.cases(data.full),]
  data.na <- data; data.na <- data.na[!complete.cases(data.na),]
  if (nrow(data.na) == 0) return(data.full)
  r <- gstat(formula=formula, locations=locations, data=data.full, nmax=7, set=list(idp=.5))
  p <- predict(r, newdata=data.na)
  data.na[,varname] <- p$var1.pred
  return(rbind(data.full, data.na))
}

theme_raster <- function(base_size=12, base_family="", ...) {
  library(grid)
  theme_minimal(base_size=base_size, base_family=base_family) %+replace%
  theme(
    panel.background=element_rect(fill="transparent", colour=NA),
    panel.grid.minor=element_blank(), 
    panel.grid.major=element_blank(),
    plot.background=element_rect(fill="transparent", colour=NA),
    #panel.background=element_blank(),
    panel.border=element_blank(),
    #panel.grid.major=element_blank(),
    #panel.grid.minor=element_blank(),
    panel.margin=unit(0, "lines"),
    axis.line=element_blank(),
    axis.text.x=element_blank(),
    axis.text.y=element_blank(),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    axis.ticks=element_blank(),
    strip.background=element_blank(),
    #plot.margin=unit(c(0,0,-1,-1), "lines"),
    plot.margin=unit(c(0,0,-.5,-.5), "lines"),
    axis.ticks.length=unit(0,"lines"),axis.ticks.margin=unit(0,"lines"),
    #plot.margin=rep(unit(0,"null"),4),panel.margin=unit(0,"null"),axis.ticks.length=unit(0,"null"),axis.ticks.margin=unit(0,"null"),
    legend.position="none",
    ...
  )
}

theme_presentation <- function(base_size=20, base_family="", ...) {
  library(grid)
  theme_minimal(base_size=base_size, base_family=base_family) %+replace%
  theme(
    ...,
    axis.line=element_line(size=1),
    legend.position="none"
  )
}

saveFigure <- function(p, filename, width=8, height=6, dimensions, ...) {
  library(ggplot2)
  if (missing(p) | missing(filename))
    stop("Argument p or filename missing.")
  filename <- file.path(context$figuresDirectory, filename)
  message("Saving plot to ", filename, "...")
  if (!missing(dimensions)) {
    aspectRatio <- dimensions[1] / dimensions[2]
    height <- width * aspectRatio
  }
  ggsave(p, filename=filename, width=width, height=height, ...)
}

addDtDist <- function(tracksDF) {
  library(data.table)
  
  if (any(!c("date", "id", "x", "y") %in% names(tracksDF)))
    stop("Need to have x, y, date and id columns in the tracks data frame.")

  tracksDF$year <- as.POSIXlt(tracksDF$date)$year + 1900
  tracksDF$yday <- as.POSIXlt(tracksDF$date)$yday
  tracksDF$burst <- paste(tracksDF$id, tracksDF$year)
  tracksDT <- data.table(tracksDF)
  
  getDxDyDtDist <- function(date,x,y) {
    n <- length(date)
    if (n == 1) return(as.numeric(NA))
    n1 <- n-1
    dx <- c(x[2:n] - x[1:n1], NA)
    dy <- c(y[2:n] - y[1:n1], NA)
    date <- as.POSIXct(date, origin="1900-01-01 00:00:00")
    dt <- c(difftime(date[2:n], date[1:n1], units="secs"), NA)
    if (any(dt[!is.na(dt)] <= 0)) warning("dt <= 0, something wrong with the data...")
    dist <- c(euclidean(cbind(x[1:n1], y[1:n1]), cbind(x[2:n], y[2:n])), NA)
    return(list(dx=dx, dy=dy, dt=dt, dist=dist))
  }
  
  tracksDT[,c("dx","dy","dt","dist") := getDxDyDtDist(date,x,y), by=burst]
  
  return(as.data.frame(tracksDT))
}

transformDeviation <- function(mean, deviation, fun, ...) {
  x <- fun(mean, ...)
  y <- mean(abs(c(x-fun(mean-deviation, ...), x-fun(mean+deviation, ...))) / qnorm(c(0.975)))
  return(y)
}

summaryStat <- function(x, rowname=NULL) {
  y <- data.frame(mean=mean(x), sd=sd(x), "0.025quant"=quantile(x, .025), "0.5quant"=quantile(x, .5), "0.975quant"=quantile(x, .975), row.names=rowname)
  colnames(y) <- c("mean", "sd", "0.025quant", "0.5quant", "0.975quant")
  return(y)
}

rangeToKappa <- function(range) sqrt(8)/range
kappaToRange <- function(kappa) sqrt(8)/kappa
sigmaToTau <- function(sigma,kappa) 1/(sqrt(4*pi)*sigma*kappa)
tauToSigma <- function(tau,kappa) 1/(sqrt(4*pi)*tau*kappa)

getSPID <- function(x) UseMethod("getSPID", x)
getSPID.SpatialLines <- function(x) sapply(x@lines, function(x) x@ID)
getSPID.SpatialLinesDataFrame <- function(x) getSPID.SpatialLines(x)
getSPID.SpatialPolygons <- function(x) sapply(x@polygons, function(x) x@ID)
getSPID.SpatialPolygonsDataFrame <- function(x) getSPID.SpatialLines(x)

getPolygonRectangle <- function(xrange, yrange, proj4string) {
  coords <- matrix(c(
    xrange[1], yrange[1],
    xrange[1], yrange[2],
    xrange[2], yrange[2],
    xrange[2], yrange[1],
    xrange[1], yrange[1]), ncol=2, byrow=T)
  r <- SpatialPolygons(list(Polygons(list(Polygon(coords)), ID="window")), proj4string=proj4string)
  return(r)
}

hasMember <- function(lst, name) if (any(name %in% names(lst))) TRUE else FALSE

colSDs <- function(x, na.rm=T) sapply(x, sd, na.rm=na.rm)

concat <- function(..., sep="") {
  args <- list(...)
  args <- args[!sapply(args, is.null)]
  paste(args, collapse=sep)
}

getMSS <- function(scenario, nSurveyRoutes, sampleInitial=F, readHabitatIntoMemory=F, isTest=FALSE, ...) {
  context <- Context(resultDataDirectory=wd.data.results, processedDataDirectory=wd.data.processed, rawDataDirectory=wd.data.raw, scratchDirectory=wd.scratch, figuresDirectory=wd.figures)
  mss <- {
    if (substr(scenario, 1, 1) == "A") MovementSimulationScenarioA(...)$setup(context=context, nSurveyRoutes=nSurveyRoutes, isTest=isTest)
    else if (substr(scenario, 1, 1) == "B") MovementSimulationScenarioB(...)$setup(context=context, nSurveyRoutes=nSurveyRoutes, isTest=isTest)
    else if (substr(scenario, 1, 1) == "C") MovementSimulationScenarioC(...)$setup(context=context, nSurveyRoutes=nSurveyRoutes, isTest=isTest)
    else if (substr(scenario, 1, 1) == "D") MovementSimulationScenarioD(...)$setup(context=context, nSurveyRoutes=nSurveyRoutes, sampleInitial=sampleInitial, isTest=isTest)
    else if (substr(scenario, 1, 1) == "E") MovementSimulationScenarioE(...)$setup(context=context, nSurveyRoutes=nSurveyRoutes, readHabitatIntoMemory=readHabitatIntoMemory, isTest=isTest)
    else if (substr(scenario, 1, 1) == "F") MovementSimulationScenarioF(...)$setup(context=context, nSurveyRoutes=nSurveyRoutes, readHabitatIntoMemory=readHabitatIntoMemory, isTest=isTest)
    else stop("unsupported")
  }
  mss$study$response <- scenario
  return(mss)
}

getStudy <- function(scenario, withHabitatWeights=FALSE, surveyRoutes, isTest=FALSE) {
  context <- Context(resultDataDirectory=wd.data.results, processedDataDirectory=wd.data.processed, rawDataDirectory=wd.data.raw, scratchDirectory=wd.scratch, figuresDirectory=wd.figures)
  study <- SimulationStudy(response=scenario)$setup(context=context, withHabitatWeights=withHabitatWeights, surveyRoutes=surveyRoutes, isTest=isTest)
  return(study)
}

parseArguments <- function() {
  args <- commandArgs(trailingOnly=TRUE)
  if (is.null(args) | length(args) == 0) {
    isTest <<- FALSE
    scenario <<- "A"
    task_id <<- as.integer(1)
    extraArgs <<- NULL
  }
  else {
    if (length(args) < 3) stop("Invalid arguments.")
    message("Arguments provided:")
    print(args)
    isTest <<- args[1] == "test"
    scenario <<- args[2]
    extraArgs <<- if (length(args) >= 3) args[3:(length(args)-1)] else NULL
    task_id <<- as.integer(args[length(args)])
  }
}
