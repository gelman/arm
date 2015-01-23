
setGeneric("coef")
setGeneric("print")
setGeneric("fitted")


#' @rdname coefplot
#' @export
setGeneric("coefplot", function(object, ...) {
  standardGeneric("coefplot")
})

#' @rdname display
#' @export
setGeneric("display", function(object, ...) {
  standardGeneric("display")
})

#' @rdname sim
#' @export
setGeneric("sim", function(object, ...) {
  standardGeneric("sim")
})

#' @rdname sigma.hat
#' @export
setGeneric("sigma.hat", function(object, ...) {
  standardGeneric("sigma.hat")
})

#' @rdname se.coef
#' @export
setGeneric("se.coef", function(object, ...) {
  standardGeneric("se.coef")
})

if (!isGeneric("fixef")) {
    setGeneric("fixef",
               function(object, ...)
               standardGeneric("fixef"),
               useAsDefault = function(object, ...) nlme::fixef(object, ...))
} 

if (!isGeneric("ranef")) {
    setGeneric("ranef",
               function(object, ...)
               standardGeneric("ranef"),
               useAsDefault = function(object, ...) nlme::ranef(object, ...))
} 

#' @rdname mcsamp
#' @export
setGeneric("mcsamp", function(object, ...) {
  standardGeneric("mcsamp")
})

#' @rdname standardize
#' @export
setGeneric("standardize", function(object, ...) {
  standardGeneric("standardize")
})

#' @rdname traceplot
#' @export
setGeneric("traceplot", function(x, ...) {
  standardGeneric("traceplot")
}, useAsDefault = function(x, ...) coda::traceplot(x, ...))


#if (!isGeneric("model.matrix.bayes")) {
#    setGeneric("model.matrix.bayes",
#               function(object, ...)
#               standardGeneric("model.matrix.bayes"))
#}
#

#if (!isGeneric("bayesglm")) {
#    setGeneric("bayesglm",
#               function(formula, ...)
#               standardGeneric("bayesglm"))
#}

#if (!isGeneric("terms.bayes")) {
#    setGeneric("terms.bayes",
#               function(x, ...)
#               standardGeneric("terms.bayes"))
#}
