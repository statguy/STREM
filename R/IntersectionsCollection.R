IntersectionsCollection <- setRefClass(
  Class = "IntersectionsCollection",
  fields = list(
    study = "Study",
    intersectionsList = "list"
  ),
  methods = list(
    initialize = function(...) {
      callSuper(...)
    },
    
    length = function() return(base::length(intersectionsList)),
    
    add = function(intersections) {
      intersectionsList[[length() + 1]] <<- intersections
    },
    
    get = function(index) {
      if (base::length(intersectionsList) < index) stop("Invalid index = ", index, ".")
      return(intersectionsList[[index]])
    },
    
    load = function() {
      stop("Unimplemented method.")
    },
    
    apply = function(fun, ...) {
      return(lapply(X=intersectionsList, FUN=fun, ...))
    }
  )
)

SimulatedIntersectionsCollection <- setRefClass(
  Class = "SimulatedIntersectionsCollection",
  contains = "IntersectionsCollection",
  methods = list(
    load = function() {     
      intersectionFiles <- context$listLongFiles(dir=study$context$resultDataDirectory, name="Intersections", response=study$response, tag="\\d+", region=study$studyArea$region)
      
      for (iteration in 1:base::length(intersectionFiles)) { # TODO: better to use iteration number rather than number of files
        intersections <- SimulatedIntersections$new(study=study, iteration=iteration)
        intersections$loadIntersections()
        add(intersections)
      }
    }
  )
)
