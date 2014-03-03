library(WTC)

# For testing

context <- Context$new(resultDataDirectory=wd.data.results, processedDataDirectory=wd.data.processed, rawDataDirectory=wd.data.raw, scratchDirectory=wd.scratch, figuresDirectory=wd.figures)
study <- FinlandWTCStudy$new(context=context, response="lynx.lynx")
model <- study$estimate(quick=TRUE)
surveyRoutes <- study$loadSurveyRoutes()
model$plotMesh(surveyRoutes)
model <- study$estimate(test=TRUE)
model <- study$estimate(quick=TRUE)

# For full estimation

context <- Context$new(resultDataDirectory=wd.data.results, processedDataDirectory=wd.data.processed, rawDataDirectory=wd.data.raw, scratchDirectory=wd.scratch, figuresDirectory=wd.figures)
study <- FinlandWTCStudy$new(context=context, response="lynx.lynx")
model <- study$estimate()
study <- FinlandWTCStudy$new(context=context, response="canis.lupus")
model <- study$estimate()
study <- FinlandWTCStudy$new(context=context, response="rangifer.tarandus.fennicus")
model <- study$estimate()
