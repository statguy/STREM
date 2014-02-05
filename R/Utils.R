euclidean <- function(c1, c2) sqrt((c1[,1] - c2[,1])^2 + (c1[,2] - c2[,2])^2)

getVector <- function(coords, distance, angle) {
  return(cbind(coords[,1] + distance * cos(angle), coords[,2] + distance * sin(angle)))
}

rotate <- function(x, y, angle) {
  return(cbind(x * cos(angle) - y * sin(angle), x * sin(angle) + y * cos(angle)))
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
