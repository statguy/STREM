

HabitatWeights <- setRefClass(
  Class = "HabitatWeights",
  fields = list(
    study = "Study",
    weightsRasterCache = "ANY",
  ),
  methods = list(
    initialize = function(...) {
      callSuper(...)
      weightsRasterCache <- environment()
      return(invisible(.self))
    },
    
    classify = function(habitatValues) {
      return(rep(1, times=length(habitatValues)))
    },
    
    getWeights = function(habitatValues) {
      return(rep(1, times=length(habitatValues)))
    },
    
    getWeightsRasterFileName = function() {
      study$context$getFileName(dir=study$context$scratchDirectory, name="HabitatWeightsRaster", response=study$response, region=study$studyArea$region, ext="")
    },
    
    getWeightsRaster = function(habitat=study$getTemplateRaster(), aggregationScale=100, save=FALSE) { # TODO: determine aggregation scale automatically
      if (save) stop("Saving unsupported.")
      library(raster)
      weightsRaster <- habitat
      weightsRaster[] <- 1
      return(weightsRaster)
    },
    
    loadWeightsRaster = function(fileName=getWeightsRasterFileName(), cache=T) {
      library(raster)
      if (cache) {
        if (exists(".weightRaster", envir=.self$weightsRasterCache))
          return(get(".weightRaster", envir=.self$weightsRasterCache))
        else {
          r <- raster(fileName)
          assign(".weightRaster", r, envir=.self$weightsRasterCache)
          return(r)
        }
      }
      else return(raster(fileName))
    },
    
    getHabitatFrequencies = function(habitatValues) {
      stop("Implement this method in a subclass.")
    },
    
    setHabitatSelectionWeights = function(relativeHabitatUsage) {
      stop("Implement this method in a subclass.")
    },
    
    show = function() {
      cat("Habitat weights:\n1\n")
      return(invisible(.self))
    }
  )
)

CORINEHabitatWeights <- setRefClass(
  Class = "CORINEHabitatWeights",
  contains = "HabitatWeights",
  fields = list(
    weights = "data.frame",
    habitatTypes = "integer",
    iteration = "integer"
  ),
  methods = list(
    initialize = function(habitatWeightsList=list(Urban=1, Agriculture=1, Forestland=1, Peatland=1, Water=1), ...) {
      callSuper(...)
      
      lastIndex <- 1+13+4+18+6+3
      
      weights <<- data.frame(
        habitat=as.integer(0:255),
        type=as.integer(c(
          rep(0,1),  # unknown
          rep(1,13), # urban
          rep(2,4),  # agriculture
          rep(3,18), # forestland
          rep(4,6),  # peatland
          rep(5,3),  # water
          rep(0,256-lastIndex)) # undefined
        ),
        weight=c(
          rep(0,1), # unknown
          rep(habitatWeightsList$Urban,13),
          rep(habitatWeightsList$Agriculture,4),
          rep(habitatWeightsList$Forestland,18),
          rep(habitatWeightsList$Peatland,6),
          rep(habitatWeightsList$Water,3),
          rep(0,256-lastIndex) # undefined
        )
      )
      
      habitatTypes <<- unique(weights$type[weights$type != 0])
      names(habitatTypes) <<- c("Urban","Agriculture","Forestland","Peatland","Water")
      
      return(invisible(.self))
    },
    
    setHabitatSelectionWeights = function(habitatSelection) {
      w <- as.list(habitatSelection$relativeUsage / habitatSelection$relativeUsage[["Forestland"]])
      initialize(habitatWeightsList=w)
      return(invisible(.self))
    },
    
    classify = function(habitatValues, na.value=NA) {
      x <- weights$type[match(habitatValues, weights$habitat)]
      x[is.na(habitatValues) | is.na(x)] <- na.value
      return(x)
    },
    
    getWeights = function(habitatValues, na.value=NA) {
      x <- weights$weight[match(habitatValues, weights$habitat)]
      x[is.na(habitatValues) | is.na(x)] <- na.value
      return(x)
    },
    
    getHabitatFrequencies = function(habitatValues) {
      mappedValues <- classify(as.vector(habitatValues))
      mappedValues <- mappedValues[!(is.na(mappedValues))]
      p <- table(factor(mappedValues,
                        levels=habitatTypes,
                        labels=c("Urban", "Agriculture", "Forestland", "Peatland", "Water")))
      return(p)
    },
    
    getWeightsRasterFileName = function() {
      if (length(iteration) == 0)
        stop("Field 'iteration' must be specified at initialization.")
      study$context$getLongFileName(dir=study$context$scratchDirectory, name="HabitatWeightsRaster", response=study$response, region=study$studyArea$region, tag=iteration, ext=".grd")
    },
    
    getWeightsRaster = function(habitat=study$studyArea$habitat, aggregationScale=100, save=FALSE, weightsRasterFileName=getWeightsRasterFileName(), grassLocalTempDir) { # TODO: determine aggregation scale automatically
      library(raster)
      
      if (file.exists(weightsRasterFileName)) {
        message("Reading habitat weights raster from ", weightsRasterFileName)
        return(.self$loadWeightsRaster(fileName=weightsRasterFileName, cache=T)) # FIXME: ad-hoc cache=T here for bootstrap estimates
      }
      
      message("Aggregating habitat weights raster...")
      
      aggfun <- function(habitatValue, na.rm) {
        x <- getWeights(habitatValue)
        x <- x[x!=0]
        mean(x, na.rm=na.rm)
      }
      
      if (!missing(grassLocalTempDir)) {
        # habitat, aggregationScale parameters ignored

        # Processes habitat raster file with GRASS GIS in HPC (but only single process in each node):
        # 1. Create environment for GRASS
        # 2. Create local temp directory to the node
        # 3. Import habitat raster to the grass mapset
        # 4. Replace habitats with weights in the raster
        # 5. Aggregate weights raster
        # 6. Export weights raster to scratch file
        # 7. Read raster to R
        # 8. Remove local GRASS environment and the scratch raster
        # There is some overhead but still faster than raster::aggregate
        
        grassCall <- paste("grass70", file.path(grassLocalTempDir, "/wtc/PERMANENT/"))
        
        # hack to parallelize grass
        grassParallel <- paste0(
"LOCALTMP=\"", grassLocalTempDir, "\"
mkdir $LOCALTMP
MYGISRC=$LOCALTMP/.grassrc
MYLOC=loc
MYMAPSET=mapset
mkdir $LOCALTMP/$MYLOC
mkdir $LOCALTMP/$MYLOC/$MYMAPSET
echo \"GISDBASE: $LOCALTMP\" > \"$MYGISRC\"
echo \"LOCATION_NAME: $MYLOC\" >> \"$MYGISRC\"
echo \"MAPSET: $MYMAPSET\" >> \"$MYGISRC\"
echo \"GRASS_GUI: text\" >> \"$MYGISRC\"
export GISRC=$MYGISRC
export TGISDB_DATABASE=$LOCALTMP/$MYLOC/PERMANENT/tgis/sqlite.db
r.in.gdal input=", study$studyArea$habitat@file@name, " output=habitat location=wtc")
        cat(grassParallel)
        system(grassParallel)
        
        
        grassRecodeInput <- ""
        for (i in 1:nrow(weights)) {
          grassRecodeInput <- paste(grassRecodeInput,
                                  paste0(weights$habitat[i],":",weights$habitat[i],":",weights$weight[i],":",weights$weight[i]),
                                  sep="\n")
        }
        grassRecodeInput <- substr(grassRecodeInput, 2, nchar(grassRecodeInput))
        
        # Note: you have must set up grass to work in text mode and the input raster must be named as 'habitat' in the mapset
        output_raster <- paste("tmp_weighted_habitat", scenario, iteration, sep="_")
        output_agg_raster <- paste("tmp_agg_weighted_habitat", scenario, iteration, sep="_")
        output_file <- file.path(study$context$scratchDirectory, paste0(output_raster, ".tif"))
        grassInput <- paste(paste0(grassCall, " << EOF"),
                            "g.region raster=habitat",
                            paste0("r.recode input=habitat output=", output_raster, " rules=-"),
                            grassRecodeInput,
                            "end",
                            "EOF", # why execution stops here?? quickfix: add second call to grass
                            paste0(grassCall, " << EOF"),
                            "g.region res=2500",
                            paste0("r.resamp.stats input=", output_raster, " output=", output_agg_raster),
                            paste0("r.out.gdal --o input=", output_agg_raster, " output=", output_file, " format=GTiff"),
                            paste0("g.remove -f type=raster name=", output_raster),
                            paste0("g.remove -f type=raster name=", output_agg_raster),
                            "exit",
                            "EOF\n",
                            sep="\n")
        cat(grassInput)
        system(grassInput)
        
        weightsRaster <- raster(output_file)
        weightsRaster <- readAll(weightsRaster)

        if (save) {
          fileName <- getWeightsRasterFileName()
          writeRaster(weightsRaster, fileName, format="raster")
        }
        
        message("Removing ", output_file)
        file.remove(output_file)
        message("Removing ", grassLocalTempDir)
        unlink(grassLocalTempDir, recursive=TRUE, force=TRUE)
      }
      else if (save) {
        fileName <- getWeightsRasterFileName()
        
        message("Saving habitat weights raster to ", fileName, ".grd...")
        weightsRaster <- aggregate(habitat,
                                   fact=aggregationScale, # TODO: get this from template raster
                                   filename=fileName, overwrite=TRUE,
                                   fun=aggfun, na.rm=T)
      }
      else {
        message("Aggregating habitat weights raster...")
        weightsRaster <- aggregate(habitat,
                                   fact=aggregationScale,
                                   fun=aggfun, na.rm=T)
      }
      
      return(weightsRaster)
    },
    
    getHabitatSelectionWeights = function() {
      w <- numeric(length(habitatTypes))
      names(w) <- names(habitatTypes)
      for (type in habitatTypes) w[type] <- weights$weight[weights$type == type][1]
      return(w)
    },
    
    show = function() {
      cat("Habitat weights:\n")
      print(getHabitatSelectionWeights())
      return(invisible(.self))
    }
  )
)
