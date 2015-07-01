WTC
===

Supplementary code for the Jousimo et al. paper
"A spatio-temporally explicit random encounter model for large-scale population surveys".

Requirements
-------------

* Linux (or possibly other UNIX alike) system HPC cluster
* 32 GB of main memory in each HPC cluster node (recommended)
* 16 CPUs in each HPC cluster (recommended)
* 500 GB of shared disk space in the HPC cluster
* 50 GB of local disk space in each HPC cluster node
* R
* Several R packages
* Python
* C/C++ compiler
* GDAL
* GEOS
* GRASS

Installation
------------

Run R and install the `devtools` package with
`install.packages("devtools")`
and the WTC R package with
`install_github("statguy/WTC")`,
FMI API package with
`install_github("statguy/WTC")` and
follow the installation instructions for the testing version R-INLA [here](http://www.r-inla.org/download).
For the rest of the dependencies, consult the installation instructions of the respective libraries.

Setup
-----

Create configuration file to your git-directory `~/git/WTC/setup/WTC-boot.R` and
add the following lines:
```
wd <- "~/wtc" # Replace with data directory
wd.figures <- file.path(wd, "figures")
wd.data <- file.path(wd, "data")
wd.data.raw <- file.path(wd.data, "raw")
wd.data.processed <- file.path(wd.data, "processed")
wd.data.results <- file.path(wd.data, "results")
wd.scratch <- file.path("/tmp/scratch") # Replace with scratch directory for large files
grassLocalTempDir <- "/tmp/wtc" # Temporary directory for GRASS
fmiApiKey <- "XYZ" # Replace with your FMI API key
```

where `wd` is the base directory for input and output data and `wd.scratch` is a scratch directory
for storing large files. `grassLocalTempDir` is a temporary directory for processing rasters with
GRASS in a HPC cluster. Replace `fmiApiKey` with your API key to the
[FMI open data](https://en.ilmatieteenlaitos.fi/open-data).

WTC uses HPC (=high-performance computing) cluster for parallel computations.
There are [a bunch of scripts](https://github.com/statguy/Parallel-R-SSH) that should be installed
from Github to use HPC. The scripts are configured for the author's HPC cluster and a custom module
for your HPC cluster should be provided.

External data installation
--------------------------

* Survey routes: copy file from `WTC/data/inst/Intersections-simulation-Finland.RData` to the directory
pointed by the variable `wd.data.processed`.
* CORINE land cover: Download the raster from http://wwwd3.ymparisto.fi/d3/Static_rs/spesific/clc2006_fi25m.zip,
unzip the file and copy the TIFF file to the directory pointed by the variable `wd.scratch` as
`HabitatRaster-Finland-cropped.tif`.
* Study area boundary polygon: Downloaded automatically from http://gadm.org.

Usage
-----

Scripts to generate simulation data, estimate models and verify estimations
are found in the directory `WTC/inst/simulations`.

Author
------

Jussi Jousimo, jvj@iki.fi
