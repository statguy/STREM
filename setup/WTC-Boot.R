wd <- "~/kolmiot"
wd.figures <- file.path(wd, "figures")
wd.data <- file.path(wd, "data")
wd.data.raw <- file.path(wd.data, "raw")
wd.data.processed <- file.path(wd.data, "processed")
wd.data.results <- file.path(wd.data, "results")
wd.scratch <- file.path("/cs/taatto/home/jousimo")

fmiApiKey <- "5bf2d697-e6df-40f5-a3d0-631e2c2642c8"

cloneData <- function(remoteMountPoint) {
  copyFiles <- function(dir, files) {
    srcDir <- if (substr(dir, 1, 1) == "~") substr(dir, 2, nchar(dir)) else dir
    src <- file.path(remoteMountPoint, srcDir)
    dest <- file.path(dir)
    
    message("Copying ", file.path(src, files), " to ", dest, "...")    
    system(paste("cp -v", file.path(src, files), dest))
  }
  
  copyFiles(wd.data.processed, "Tracks-*")
  copyFiles(wd.data.processed, "FinlandWTCIntersectionsCovariates-*")
  copyFiles(wd.data.results, "SmoothModel-*")
  copyFiles(wd.data.results, "Intersections-*")
  copyFiles(wd.data.results, "HabitatWeights-*")
  copyFiles(wd.data.results, "HabitatWeightsRaster-*")

  # Does not copy scratch data
}

#cloneData("/mnt/melkinpaasi")
