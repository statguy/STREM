PopulationSize <- setRefClass(
  Class = "PopulationSize",
  fields = list(
    size = "data.frame"
  ),
  methods = list(
    add = function(year, size) {
      size <<- rbind(size, data.frame(Year=year, Estimated=size))
    }
    
  )
)