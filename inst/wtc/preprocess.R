library(devtools)
install_github("R-Cluster", user="statguy")
install_github("R-Spatio-Temporal", user="statguy")
install_github("Winter-Track-Counts", user="statguy")

fmiApiKey <- "" # Specify your FMI API key here

library(WTC)

context <- Context$new(resultDataDirectory=wd.data.results, processedDataDirectory=wd.data.processed, rawDataDirectory=wd.data.raw, scratchDirectory=wd.scratch, figuresDirectory=wd.figures)
study <- FinlandWTCStudy$new(context=context)
study$preprocess(fmiApiKey=fmiApiKey)
