study <- getStudy("A")
transects <- FinlandRandomForestWTCSurveyRoutes$new(study=study)
transects$randomizeSurveyRoutes(500)
transects$saveSurveyRoutes()

transects$getSurveyRoutesFileName()
