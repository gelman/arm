setMethod("sigma.hat", signature(object = "lm"),
    function(object)
    {
    object.class <- class(object)[[1]]
    sigma <- summary(object)$sigma
    return (sigma)
    }
)


setMethod("sigma.hat", signature(object = "glm"),
    function(object)
    {
    dispersion <- if (is.null(object$dispersion)){
                    summary(object)$dispersion
                  }
                  else{
                    object$dispersion
                  }
    object.class <- class(object)[[1]]
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



setMethod("sigma.hat", signature(object = "mer"),
    function(object)
    {
    #object <- summary (object)
    fcoef <- fixef(object)
    #useScale <- attr (VarCorr (object), "sc")  # =sc?
    useScale <- object@dims["useSc"]
    ngrps <- lapply(object@flist, function(x) length(levels(x)))
    n.groupings <- length (ngrps)
    varc <- VarCorr (object)
    sc <- attr(varc, "sc")  # =useScale
    recorr <- lapply(varc, function(el) attr(el, "correlation"))
    reStdDev <- c(lapply(varc, function(el) attr(el, "stddev")), list(Residual = sc))
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
