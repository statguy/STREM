library(WTC)

context <- Context$new(resultDataDirectory=wd.data.results, processedDataDirectory=wd.data.processed, rawDataDirectory=wd.data.raw, scratchDirectory=wd.scratch, figuresDirectory=wd.figures)
study <- FinlandWTCStudy$new(context=context, response="lynx.lynx")
model <- study$estimate(test=TRUE)
surveyRoutes <- study$loadSurveyRoutes()
model$plotMesh(surveyRoutes)

model <- study$estimate()
