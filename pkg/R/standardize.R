#' Function for Standardizing Regression Predictors by Centering and
#' Dividing by 2 sd's
#'
#' Numeric variables that take on more than two values are each
#' rescaled to have a mean of 0 and a sd of 0.5; Binary variables are
#' rescaled to have a mean of 0 and a difference of 1 between their
#' two categories; Non-numeric variables that take on more than two
#' values are unchanged; Variables that take on only one value are
#' unchanged.
#'
#' "0/1" (rescale so that the lower value is 0 and the upper is 1)
#' "-0.5/0.5" (rescale so that the lower value is -0.5 and upper is 0.5)
#' "center" (rescale so that the mean of the data is 0 and the difference
#'           between the two categories is 1)
#' "full" (rescale by subtracting the mean and dividing by 2 sd's)
#' "leave.alone" (do nothing)
#'
#' @param object an object of class \code{lm} or \code{glm}.
#' @param unchanged vector of names of parameters to leave
#'   unstandardized.
#' @param standardize.y if TRUE, the outcome variable is standardized
#'   also.
#' @param binary.inputs options for standardizing binary variables.
#' @references Andrew Gelman. (2008). \dQuote{Scaling regression
#'   inputs by dividing by two standard deviations.} \emph{Statistics in
#'   Medicine} 27: 2865--2873.
#'   \url{http://www.stat.columbia.edu/~gelman/research/published/standardizing7.pdf}
#' @author
#'   Andrew Gelman \email{gelman@@stat.columbia.edu}
#'   Yu-Sung Su \email{suyusung@@tsinghua.edu.cn}
#' @seealso \code{\link{rescale}}
#' @keywords manip models methods
#' @name standardize
#' @examples
#' # Set up the fake data
#' n <- 100
#' x <- rnorm (n, 2, 1)
#' x1 <- rnorm (n)
#' x1 <- (x1-mean(x1))/(2*sd(x1))   # standardization
#' x2 <- rbinom (n, 1, .5)
#' b0 <- 1
#' b1 <- 1.5
#' b2 <- 2
#' y <- rbinom (n, 1, invlogit(b0+b1*x1+b2*x2))
#' y2 <- sample(1:5, n, replace=TRUE)
#' M1 <- glm (y ~ x, family=binomial(link="logit"))
#' display(M1)
#' M1.1 <- glm (y ~ rescale(x), family=binomial(link="logit"))
#' display(M1.1)
#' M1.2 <- standardize(M1.1)
#' display(M1.2)
#' # M1.1 & M1.2 should be the same
#' M2 <- polr(ordered(y2) ~ x)
#' display(M2)
#' M2.1 <- polr(ordered(y2) ~ rescale(x))
#' display(M2.1)
#' M2.2 <- standardize(M2.1)
#' display(M2.2)
#' # M2.1 & M2.2 should be the same
NULL

standardize.default <- function(call, unchanged=NULL,
                         standardize.y=FALSE, binary.inputs="center"){
  form <- call$formula
  varnames <- all.vars (form)
  n.vars <- length (varnames)
#
# Decide which variables will be unchanged
#
  transform <- rep ("leave.alone", n.vars)
  if (standardize.y) {
    transform[1] <- "full"
  }
  for (i in 2:n.vars){
    v <- varnames[i]
    if (is.null(call$data)) {
      thedata <- get(v)
    }
    else {
      thedata <- get(as.character(call$data))[[v]]
    }
    if (is.na(match(v,unchanged))){
      num.categories <- length (unique(thedata[!is.na(thedata)]))
      if (num.categories==2){
        transform[i] <- binary.inputs
      }
      else if (num.categories>2 & is.numeric(thedata)){
        transform[i] <- "full"
      }
    }
  }
#
# New variable names:
#   prefix with "c." if centered or "z." if centered and scaled
#
  varnames.new <- ifelse (transform=="leave.alone", varnames,
    ifelse (transform=="full", paste ("z", varnames, sep="."),
    paste ("c", varnames, sep=".")))
  transformed.variables <- (1:n.vars)[transform!="leave.alone"]


  #Define the new variables
  if (is.null(call$data)) {
    for (i in transformed.variables) {
      assign(varnames.new[i], rescale(get(varnames[i]), binary.inputs))
    }
  }
  else {
    newvars <- NULL
    for (i in transformed.variables) {
      assign(varnames.new[i], rescale(get(as.character(call$data))[[varnames[i]]], 
                binary.inputs))
      newvars <- cbind(newvars, get(varnames.new[i]))
    }
    assign(as.character(call$data), cbind(get(as.character(call$data)), newvars))
  }

# Now call the regression with the new variables

  call.new <- call
  L <- sapply (as.list (varnames.new), as.name)
  names(L) <- varnames
  call.new$formula <- do.call (substitute, list (form, L))
  formula <- as.character (call.new$formula)
  if (length(formula)!=3) stop ("formula does not have three components")
  formula <- paste (formula[2],formula[1],formula[3])
  formula <- gsub ("factor(z.", "factor(", formula, fixed=TRUE)
  formula <- gsub ("factor(c.", "factor(", formula, fixed=TRUE)
  call.new$formula <- as.formula (formula) 
  return (eval (call.new))
}


#' @rdname standardize
#' @export
setMethod("standardize", signature(object = "lm"),
  function(object, unchanged=NULL, 
    standardize.y=FALSE, binary.inputs="center")
{
  call <- object$call
  out <- standardize.default(call=call, unchanged=unchanged, 
    standardize.y=standardize.y, binary.inputs=binary.inputs)
  return(out)
}
)

#' @rdname standardize
#' @export
setMethod("standardize", signature(object = "glm"),
  function(object, unchanged=NULL, 
    standardize.y=FALSE, binary.inputs="center")
{
  call <- object$call
  out <- standardize.default(call=call, unchanged=unchanged, 
    standardize.y=standardize.y, binary.inputs=binary.inputs)
  return(out)
}
)

#' @rdname standardize
#' @export
setMethod("standardize", signature(object = "polr"),
  function(object, unchanged=NULL, 
    standardize.y=FALSE, binary.inputs="center")
{
  call <- object$call
  out <- standardize.default(call=call, unchanged=unchanged, 
    standardize.y=standardize.y, binary.inputs=binary.inputs)
  return(out)
}
)

#' @rdname standardize
#' @export
setMethod("standardize", signature(object = "merMod"),
  function(object, unchanged=NULL, 
    standardize.y=FALSE, binary.inputs="center")
{
  call <- object@call
  out <- standardize.default(call=call, unchanged=unchanged, 
    standardize.y=standardize.y, binary.inputs=binary.inputs)
  return(out)
}
)
