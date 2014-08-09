#if (!isGeneric("bayesglm")) {
#    setGeneric("bayesglm",
#               function(formula, ...)
#               standardGeneric("bayesglm"))
#}


setGeneric("coef")
setGeneric("print")
setGeneric("fitted")

#setGeneric("extractAIC")

if (!isGeneric("coefplot")) {
    setGeneric("coefplot",
               function(object, ...)
               standardGeneric("coefplot"))
}


if (!isGeneric("display")) {
    setGeneric("display",
               function(object, ...)
               standardGeneric("display"))
}

#if (!isGeneric("model.matrix.bayes")) {
#    setGeneric("model.matrix.bayes",
#               function(object, ...)
#               standardGeneric("model.matrix.bayes"))
#}
#

if (!isGeneric("sim")) {
    setGeneric("sim",
               function(object, ...)
               standardGeneric("sim"))
}


if (!isGeneric("sigma.hat")) {
    setGeneric("sigma.hat",
               function(object, ...)
               standardGeneric("sigma.hat"))
}

if (!isGeneric("se.coef")) {
    setGeneric("se.coef",
               function(object, ...)
               standardGeneric("se.coef"))
}

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


if (!isGeneric("mcsamp")) {
    setGeneric("mcsamp",
               function(object, ...)
               standardGeneric("mcsamp"))
}



if (!isGeneric("standardize")) {
    setGeneric("standardize",
               function(object, ...)
               standardGeneric("standardize"))
}



#if (!isGeneric("terms.bayes")) {
#    setGeneric("terms.bayes",
#               function(x, ...)
#               standardGeneric("terms.bayes"))
#}


if (!isGeneric("traceplot")) {
    setGeneric("traceplot",
               function(x, ...)
               standardGeneric("traceplot"),
               useAsDefault = function(x, ...) coda::traceplot(x, ...))
} 
