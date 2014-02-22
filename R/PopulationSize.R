PopulationSize <- setRefClass(
  Class = "PopulationSize",
  fields = list(
    sizeData = "data.frame"
  ),
  methods = list(
    initialize = function(...) {
      callSuper(sizeData=data.frame(), ...)
      return(invisible(.self))
    },
    
    addYearSize = function(year, size) {
      sizeData <<- rbind(sizeData, data.frame(Year=year, Estimated=size))
    }
    
  )
)
