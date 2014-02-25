library(WTC)

fmiApiKey <- "" # Specify your FMI API key here

context <- Context$new(resultDataDirectory=wd.data.results, processedDataDirectory=wd.data.processed, rawDataDirectory=wd.data.raw, scratchDirectory=wd.scratch, figuresDirectory=wd.figures)
study <- FinlandWTCStudy$new(context=context)
study$preprocess(fmiApiKey)
