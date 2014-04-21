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
    plot.margin=unit(c(0,0,-1,-1), "lines"),
    legend.position="none",
    ...
  )
}

theme_presentation <- function(base_size=20, base_family="", ...) {
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

NEW_addDtDist <- function(tracksDF, .parallel=FALSE) {
  library(dplyr)
  
  if (any(!c("date", "id", "x", "y") %in% names(tracksDF)))
    stop("Need to have x, y, date and id columns in the tracks data frame.")
  tracksDF$year <- as.POSIXlt(tracksDF$date)$year + 1900
  tracksDF$yday <- as.POSIXlt(tracksDF$date)$yday
  tracksDF$burst <- paste(tracksDF$id, tracksDF$year, tracksDF$yday)
  
  addDt <- function(date) {
    n <- length(date)
    if (n == 1) return(as.numeric(NA))
    date <- as.POSIXct(date, origin="1900-01-01 00:00:00")
    n1 <- n-1
    dt <- c(difftime(date[2:n], date[1:n1], units="secs"), NA)
    if (any(dt[!is.na(dt)] <= 0)) warning("dt <= 0, something wrong with the data...")
    return(dt)
  }
  
  addDist <- function(x, y) {
    n <- length(x)
    if (n == 1) return(as.numeric(NA))
    n1 <- n-1
    dist <- c(euclidean(cbind(x[1:n1], y[1:n1]), cbind(x[2:n], y[2:n])), NA)
    return(dist)
  }
  
  x <- tracksDF %.%
    group_by(id, year, yday) %.%
    mutate(dt=addDt(date), dist=addDist(x,y))
  
  return(x)
}

addDtDist <- function(tracksDF, .parallel=TRUE) {
  library(data.table)
  
  if (any(!c("date", "id", "x", "y") %in% names(tracksDF)))
    stop("Need to have x, y, date and id columns in the tracks data frame.")

  tracksDF$year <- as.POSIXlt(tracksDF$date)$year + 1900
  tracksDF$yday <- as.POSIXlt(tracksDF$date)$yday
  tracksDF$burst <- paste(tracksDF$id, tracksDF$year, tracksDF$yday)
  tracksDT <- data.table(tracksDF)
  
  getDtDist <- function(date,x,y) {
    n <- length(date)
    if (n == 1) return(as.numeric(NA))
    date <- as.POSIXct(date, origin="1900-01-01 00:00:00")
    n1 <- n-1
    dt <- c(difftime(date[2:n], date[1:n1], units="secs"), NA)
    if (any(dt[!is.na(dt)] <= 0)) warning("dt <= 0, something wrong with the data...")
    dist <- c(euclidean(cbind(x[1:n1], y[1:n1]), cbind(x[2:n], y[2:n])), NA)
    return(list(dt=dt, dist=dist))
  }
  
  tracksDT[,c("dt","dist") := getDtDist(date,x,y), by=burst]
  
  return(as.data.frame(tracksDT))
}

OLD_addDtDist <- function(tracksDF, .parallel=TRUE) {
  library(plyr)
  
  if (any(!c("date", "id", "x", "y") %in% names(tracksDF)))
    stop("Need to have x, y, date and id columns in the tracks data frame.")
  tracksDF$year <- as.POSIXlt(tracksDF$date)$year + 1900
  tracksDF$yday <- as.POSIXlt(tracksDF$date)$yday
  tracksDF$burst <- paste(tracksDF$id, tracksDF$year, tracksDF$yday)
  
  tracksDF <- ddply(tracksDF, .(burst), function(x) {
    n <- nrow(x)
    if (n == 1) {
      x$dt <- NA
      x$dist <- NA
      return(x)
    }
    n1 <- n-1
    x$dt <- c(difftime(x$date[2:n], x$date[1:n1], units="secs"), NA)
    x$dist <- c(euclidean(x[1:n1, c("x","y")], x[2:n, c("x","y")]), NA)
    if (any(x$dt[!is.na(x$dt)] <= 0)) warning("dt <= 0, something wrong with the data...")
    return(x)
  }, .parallel=.parallel, .progress="text")
  
  return(tracksDF)
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
