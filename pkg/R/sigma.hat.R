#' Extract Residual Errors
#'
#' This generic function extracts residual errors from a fitted model.
#'
#' @param object any fitted model object of \code{lm}, \code{glm} and
#'   \code{merMod} class.
#' @param ... other arguments.
#' @author
#'   Andrew Gelman \email{gelman@@stat.columbia.edu};
#'   Yu-Sung Su \email{suyusung@@tsinghua.edu.cn}
#' @seealso \code{\link{display}}, \code{\link{summary}},
#' \code{\link{lm}}, \code{\link{glm}}, \code{\link[lme4]{lmer}}
#' @keywords manip methods
#' @name sigma.hat
#' @examples
#' library("MASS")
#' library("lme4")
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
#' M1 <- lm (y1 ~ x) 
#' sigma.hat(M1)
#' 
#' M2 <- bayesglm (y1 ~ x, prior.scale=Inf, prior.df=Inf)
#' sigma.hat(M2) # should be same to sigma.hat(M1)
#'
#' M3 <- glm (y2 ~ x, family=binomial(link="logit"))
#' sigma.hat(M3)
#'
#' M4 <- lmer (y1 ~ (1+x|group))
#' sigma.hat(M4) 
#' 
#' M5 <- glmer (y2 ~ (1+x|group), family=binomial(link="logit"))
#' sigma.hat(M5)
NULL

#' @rdname sigma.hat
#' @export
setMethod("sigma.hat", signature(object = "lm"),
    function(object)
    {
    sigma <- summary(object)$sigma
    return (sigma)
    }
)

#' @rdname sigma.hat
#' @export
setMethod("sigma.hat", signature(object = "sim"),
    function(object)
    {
    sigma <- object@sigma
    return (sigma)
    }
)

#' @rdname sigma.hat
#' @export
setMethod("sigma.hat", signature(object = "glm"),
    function(object)
    {
    dispersion <- if (is.null(object$dispersion)){
                    summary(object)$dispersion
                  }
                  else{
                    object$dispersion
                  }
    if (object$family$family == "gaussian") {
      sigma <- sqrt(dispersion)
    }
    else {
      sigma <- summary(object, correlation = TRUE)$sigma
      #sigma <- sqrt(deviance(object)/df.residual(object))
    }
    return(sigma)
    }
)

#setMethod("sigma.hat", signature(object = "mer"),
#    function(object)
#    {
#    #object <- summary (object)
#    fcoef <- fixef(object)
#    useScale <- object@devComp[8]
#    ngrps <- lapply(object@flist, function(x) length(levels(x)))
#    n.groupings <- length (ngrps)
#    varc <- VarCorr (object, useScale=useScale)
#    sc <- attr(varc, "sc")
#    recorr <- lapply(varc, function(el) el@factors$correlation)
#    reStdDev <- c(lapply(recorr, slot, "sd"), list(Residual = sc))
#    sigmas <- as.list (rep (NA, n.groupings+1))
#    sigmas[1] <- ifelse (useScale, sc, NA)
#    cors <- as.list (rep (NA, n.groupings+1))
#    names (sigmas) <- names (cors) <- c ("data", names (varc))
#    for (k in 1:n.groupings){
#      sigmas[[k+1]] <- reStdDev[[k]]
#      cors[[k+1]] <- as.matrix (recorr[[k]])
#      if (length (cors[[k+1]]) == 1) cors[[k+1]] <- NA
#    }
#    return (list (sigma=sigmas, cors=cors))
#    }
#) 

#' @rdname sigma.hat
#' @export
setMethod("sigma.hat", signature(object = "merMod"),
    function(object)
    {
    #object <- summary (object)
    fcoef <- fixef(object)
    #useScale <- attr (VarCorr (object), "sc")  # =sc?
    #useScale <- object@dims["useSc"]
    useScale <- getME(object, "devcomp")$dims["useSc"]
    #ngrps <- lapply(object@flist, function(x) length(levels(x)))
    #n.groupings <- length (ngrps)
    varc <- VarCorr (object)
    sc <- attr(varc, "sc")  # =useScale
    recorr <- lapply(varc, function(el) attr(el, "correlation"))
    reStdDev <- c(lapply(varc, function(el) attr(el, "stddev")), list(Residual = sc))
    n.groupings <- length(recorr)
    sigmas <- as.list (rep (NA, n.groupings+1))
    sigmas[1] <- ifelse (useScale, sc, 1) #####if NA, sd=1
    cors <- as.list (rep (NA, n.groupings+1))
    names (sigmas) <- names (cors) <- c ("data", names (varc))
    for (k in 1:n.groupings){
      sigmas[[k+1]] <- reStdDev[[k]]
      cors[[k+1]] <- as.matrix (recorr[[k]])
      if (length (cors[[k+1]]) == 1) cors[[k+1]] <- NA
    }
    return (list (sigma=sigmas, cors=cors))
    }
)  

#' @rdname sigma.hat
#' @export
setMethod("sigma.hat", signature(object = "sim.merMod"),
    function(object)
    {
    sigma <- object@sigma
    return (sigma)
    }
)
