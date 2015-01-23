#' Extract Standard Errors of Model Coefficients
#'
#' These functions extract standard errors of model coefficients from
#' objects returned by modeling functions.
#'
#' \code{se.coef} extracts standard errors from objects returned by
#' modeling functions. \code{se.fixef} extracts standard errors of the
#' fixed effects from objects returned by lmer and glmer
#' functions. \code{se.ranef} extracts standard errors of the random
#' effects from objects returned by lmer and glmer functions.
#' 
#' @param object object of \code{lm}, \code{glm} and \code{merMod}
#' fit.
#' @param ... other arguments.
#' @return \code{se.coef} gives lists of standard errors for
#' \code{coef}, \code{se.fixef} gives a vector of standard errors for
#' \code{fixef} and \code{se.ranef} gives a list of standard errors
#' for \code{ranef}.
#' @references Andrew Gelman and Jennifer Hill. (2006). \emph{Data
#' Analysis Using Regression and Multilevel/Hierarchical
#' Models}. Cambridge University Press.
#' @author
#'   Andrew Gelman \email{gelman@@stat.columbia.edu};
#'   Yu-Sung Su \email{suyusung@@tsinghua.edu.cn}
#' @seealso \code{\link{display}}, \code{\link{coef}},
#' \code{\link{sigma.hat}}
#' @keywords manip methods models
#' @name se.coef
#' @examples
#' # Here's a simple example of a model of the form, y = a + bx + error, 
#' # with 10 observations in each of 10 groups, and with both the 
#' # intercept and the slope varying by group.  First we set up the model and data.
#'
#' group <- rep(1:10, rep(10,10))
#' mu.a <- 0
#' sigma.a <- 2
#' mu.b <- 3
#' sigma.b <- 4
#' rho <- 0
#' Sigma.ab <- array (c(sigma.a^2, rho*sigma.a*sigma.b, 
#'                  rho*sigma.a*sigma.b, sigma.b^2), c(2,2))
#' sigma.y <- 1
#' ab <- mvrnorm (10, c(mu.a,mu.b), Sigma.ab)
#' a <- ab[,1]
#' b <- ab[,2]
#'
#' x <- rnorm (100)
#' y1 <- rnorm (100, a[group] + b[group]*x, sigma.y)
#' y2 <- rbinom(100, 1, prob=invlogit(a[group] + b*x))
#'
#' # lm fit
#' M1 <- lm (y1 ~ x)
#' se.coef (M1)
#' 
#' # glm fit
#' M2 <- glm (y2 ~ x)
#' se.coef (M2)
#'
#' # lmer fit
#' M3 <- lmer (y1 ~ x + (1 + x |group))
#' se.coef (M3)
#' se.fixef (M3)
#' se.ranef (M3)
#'
#' # glmer fit
#' M4 <- glmer (y2 ~ 1 + (0 + x |group), family=binomial(link="logit"))
#' se.coef (M4)
#' se.fixef (M4)
#' se.ranef (M4)
NULL

#' @rdname se.coef
#' @export
setMethod("se.coef", signature(object = "lm"),
    function(object)
    {
    object.class <- class(object)[[1]]
    sqrt (diag(vcov(object)))
    }
)

#' @rdname se.coef
#' @export
setMethod("se.coef", signature(object = "glm"),
    function(object)
    {
    object.class <- class(object)[[1]]
    sqrt (diag(vcov(object)))
    }
)

#setMethod("se.coef", signature(object = "mer"),
#    function(object)
#    {
#    #    if (sum(unlist(lapply(object@bVar, is.na)))>0){
##        object@call$control <- list(usePQL=TRUE)
##        object <- lmer(object@call$formula)
##    }
#    #ngrps <- lapply(object@flist, function(x) length(levels(x)))
#    fcoef <- fixef(object)
#    #sc <- attr (VarCorr (object), "sc")
#    corF <- vcov(object)@factors$correlation
#    se.unmodeled <- NULL
#    se.unmodeled[[1]] <- corF@sd
#    names (se.unmodeled) <- "unmodeled"
#
#    #coef <- ranef (object)
#    #estimate <- ranef(object, postVar=TRUE)
#    coef <- ranef(object, postVar=TRUE)
#    se.bygroup <- coef #ranef( object, postVar = TRUE )
#    n.groupings <- length (coef)
#    
#    for (m in 1:n.groupings){
#      vars.m <- attr (coef[[m]], "postVar")
#      K <- dim(vars.m)[1]
#      J <- dim(vars.m)[3]
#      se.bygroup[[m]] <- array (NA, c(J,K))
#      for (j in 1:J){
#        se.bygroup[[m]][j,] <- sqrt(diag(as.matrix(vars.m[,,j])))
#      }
##      se.bygroup[[m]] <- se.bygroup[[m]]*sc
#      names.full <- dimnames (ranef(object)[[m]])
#      dimnames (se.bygroup[[m]]) <- list (names.full[[1]],
#                            names.full[[2]])
#    }
#    #names(se.bygroup) <- names(ngrps)
#    ses <- c (se.unmodeled, se.bygroup)
#    return (ses)
#    }
#)

#' @rdname se.coef
#' @export
setMethod("se.coef", signature(object = "merMod"),
    function(object)
    {
    #ngrps <- lapply(object@flist, function(x) length(levels(x)))
    fcoef <- fixef(object)
    #sc <- attr (VarCorr (object), "sc")
    corF <- vcov(object)@factors$correlation
    se.unmodeled <- NULL
    se.unmodeled[[1]] <- corF@sd
    names (se.unmodeled) <- "fixef"#"unmodeled"

    #coef <- ranef (object)
    #estimate <- ranef(object, postVar=TRUE)
    coef <- ranef(object, condVar=TRUE)
    se.bygroup <- coef #ranef( object, postVar = TRUE )
    n.groupings <- length (coef)
    
    for (m in 1:n.groupings){
      vars.m <- attr (coef[[m]], "postVar")
      K <- dim(vars.m)[1]
      J <- dim(vars.m)[3]
      se.bygroup[[m]] <- array (NA, c(J,K))
      for (j in 1:J){
        se.bygroup[[m]][j,] <- sqrt(diag(as.matrix(vars.m[,,j])))
      }
#      se.bygroup[[m]] <- se.bygroup[[m]]*sc
      names.full <- dimnames (coef[[m]])
      dimnames (se.bygroup[[m]]) <- list (names.full[[1]],
                            names.full[[2]])
    }
    #names(se.bygroup) <- names(ngrps)
    ses <- c (se.unmodeled, se.bygroup)
    return (ses)
    }
)

#' @rdname se.coef
#' @export
se.fixef <- function (object){
  #object <- summary (object)
  fcoef.name <- names(fixef(object))
  corF <- vcov(object)@factors$correlation
  ses <- corF@sd
  names(ses) <- fcoef.name
  return (ses)
}

#' @rdname se.coef
#' @export
se.ranef <- function (object){
    #ngrps <- lapply(object@flist, function(x) length(levels(x)))
    se.bygroup <- ranef( object, condVar = TRUE )
    n.groupings<- length( se.bygroup )
    for( m in 1:n.groupings ) {
        vars.m <- attr( se.bygroup[[m]], "postVar" )
        K <- dim(vars.m)[1]
        J <- dim(vars.m)[3]
        names.full <- dimnames(se.bygroup[[m]])
        se.bygroup[[m]] <- array(NA, c(J, K))
        for (j in 1:J) {
            se.bygroup[[m]][j, ] <- sqrt(diag(as.matrix(vars.m[, , j])))
        }       
        dimnames(se.bygroup[[m]]) <- list(names.full[[1]], names.full[[2]])
    }
    return(se.bygroup)
}
