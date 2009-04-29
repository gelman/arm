setMethod("display", signature(object = "lm"),
    function(object, digits=2, detail=FALSE)
    {
    call <- object$call
    summ <- summary (object)
    if(detail){
      coef <- summ$coef[,,drop=FALSE]
    }
    else{
      coef <- summ$coef[,1:2,drop=FALSE]
    }
    dimnames(coef)[[2]][1:2] <- c("coef.est","coef.se")
    n <- summ$df[1] + summ$df[2]
    k <- summ$df[1]
    print (call)
    pfround (coef, digits)
    cat("---\n")
    cat (paste ("n = ", n, ", k = ", k,
    "\nresidual sd = ", fround (summ$sigma, digits),
    ", R-Squared = ", fround (summ$r.squared, 2), "\n", sep=""))
    }
)



setMethod("display", signature(object = "bayesglm"),
    function(object, digits=2, detail=FALSE)
    {
    call <- object$call
    summ <- summary(object, dispersion = object$dispersion)
    if(detail){
      coef <- summ$coefficients
      coef[ rownames( coef ) %in% rownames( summ$coef[, , drop = FALSE]) , ] <- summ$coef[ , , drop = FALSE ] 
    }
    else{
      coef <- matrix( NA, length( object$coefficients ),2 )
      rownames(coef) <- names( object$coefficients )          ## M
      coef[ rownames( coef ) %in% rownames( summ$coef[, 1:2, drop = FALSE]) , ] <- summ$coef[ , 1:2, drop = FALSE ]  ## M
    }
    dimnames(coef)[[2]][1:2] <- c( "coef.est", "coef.se")
    #n <- summ$df[1] + summ$df[2]
    n <- summ$df.residual
    k <- summ$df[1]
    print(call)
    pfround(coef, digits)
    cat("---\n")
    cat(paste("n = ", n, ", k = ", k, "\nresidual deviance = ", 
        fround(summ$deviance, 1), ", null deviance = ", fround(summ$null.deviance, 
            1), " (difference = ", fround(summ$null.deviance - 
            summ$deviance, 1), ")", "\n", sep = ""))
    dispersion <- if (is.null(object$dispersion)) 
        summ$dispersion
    else object$dispersion
    if (dispersion != 1) {
        cat(paste("overdispersion parameter = ", fround(dispersion, 
            1), "\n", sep = ""))
        if (family(object)$family == "gaussian") {
            cat(paste("residual sd is sqrt(overdispersion) = ", 
                fround(sqrt(dispersion), digits), "\n", sep = ""))
        }
    }
    }
)

setMethod("display", signature(object = "bayesglm.h"),
    function (object, digits = 2) 
    {
    call <- object$call
    summ <- summary(object, dispersion = object$dispersion)
    if(detail){
      coef <- summ$coefficients
      coef[ rownames( coef ) %in% rownames( summ$coef[, , drop = FALSE]) , ] <- summ$coef[ , , drop = FALSE ] 
    }
    else{
      coef <- matrix( NA, length( object$coefficients ),2 )
      rownames(coef) <- names( object$coefficients )          ## M
      coef[ rownames( coef ) %in% rownames( summ$coef[, 1:2, drop = FALSE]) , ] <- summ$coef[ , 1:2, drop = FALSE ]  ## M
    }
    dimnames(coef)[[2]][1:2] <- c( "coef.est", "coef.se")
    #n <- summ$df[1] + summ$df[2]
    n <- summ$df.residual
    k <- summ$df[1]
    print(call)
    if(max(object$batch)>0){
        nn<- strsplit( rownames( coef )[seq( from= length( object$batch ) + 1 ,to = nrow( coef ))], "." , fixed=TRUE)
        bb<- c( object$batch,unlist( lapply (nn , function( lst ) { lst[[3]] } ) ) ) 
    }
    else {bb<- c( object$batch)}
    cc<- cbind( fround( coef, digits ), bb )
    dimnames(cc)[[2]][3]<-"batch"
    print( cc , quote = FALSE )
    cat("---\n")
    cat(paste("n = ", n, ", k = ", k, "\nresidual deviance = ", 
        fround(summ$deviance, 1), ", null deviance = ", fround(summ$null.deviance, 
            1), " (difference = ", fround(summ$null.deviance - 
            summ$deviance, 1), ")", "\n", sep = ""))
    dispersion <- if (is.null(object$dispersion)) 
        summ$dispersion
    else object$dispersion
    if (dispersion != 1) {
        cat(paste("overdispersion parameter = ", fround(dispersion, 
            1), "\n", sep = ""))
        if (family(object)$family == "gaussian") {
            cat(paste("residual sd is sqrt(overdispersion) = ", 
                fround(sqrt(dispersion), digits), "\n", sep = ""))
            cat(paste("group sd is sigma.batch = ", 
                fround(object$sigma.batch, digits), "\n", sep = ""))
        }
    }
    }
)




setMethod("display", signature(object = "glm"),
    function(object, digits=2, detail=FALSE)
    {
    call <- object$call
    summ <- summary(object, dispersion = object$dispersion)
    if(detail){
      coef <- summ$coef[, , drop = FALSE]
    }
    else{
      coef <- summ$coef[, 1:2, drop = FALSE]
    }
    dimnames(coef)[[2]][1:2] <- c("coef.est", "coef.se")
    n <- summ$df[1] + summ$df[2]
    k <- summ$df[1]
    print(call)
    pfround(coef, digits)
    cat("---\n")
    cat(paste("  n = ", n, ", k = ", k, "\n  residual deviance = ", 
        fround(summ$deviance, 1), ", null deviance = ", fround(summ$null.deviance, 
            1), " (difference = ", fround(summ$null.deviance - 
            summ$deviance, 1), ")", "\n", sep = ""))
    dispersion <- if (is.null(object$dispersion))
      summ$dispersion
    else object$dispersion
    if (dispersion != 1) {
      cat(paste("  overdispersion parameter = ",
                fround(dispersion, 1), "\n", sep = ""))
      if (family(object)$family=="gaussian") {
        cat(paste("  residual sd is sqrt(overdispersion) = ",
                  fround(sqrt(dispersion), digits), "\n", sep = ""))
      }
    }
    }
)




#setMethod("display", signature(object = "mer"),
#    function(object, digits=2)
#    {
#    call <- object@call
#    print (call)
#    #object <- summary(object)
#    fcoef <- fixef(object)
#    useScale <- attr( VarCorr(object), "sc") 
#    corF <- vcov(object)@factors$correlation
#    coefs <- cbind(fcoef, corF@sd)
#    if (length (fcoef) > 0){
#      dimnames(coefs) <- list(names(fcoef), c("coef.est", "coef.se"))
#      pfround (coefs, digits)
#    }
#    cat("\nError terms:\n")
#    vc <- as.matrix.VarCorr (VarCorr (object), useScale=useScale, digits)
#    print (vc[,c(1:2,4:ncol(vc))], quote=FALSE)
#    ngrps <- lapply(object@flist, function(x) length(levels(x)))
#    REML <- object@status["REML"]
#    llik <- logLik(object, REML)
#    AIC <- AIC(llik)
#    dev <- object@deviance["ML"]     # Dbar
#    n <- object@devComp["n"]
#    Dhat <- -2*(llik) # Dhat
#    pD <- dev - Dhat              # pD
#    DIC <- dev + pD               # DIC=Dbar+pD=Dhat+2pD
#    cat("---\n")
#    cat(sprintf("number of obs: %d, groups: ", n))
#    cat(paste(paste(names(ngrps), ngrps, sep = ", "), collapse = "; "))
#    cat(sprintf("\nAIC = %g, DIC = ", fround(AIC, 1)))
#    cat(fround(DIC, 1))
#    cat("\ndeviance =", fround (dev, 1), "\n")
#    if (useScale < 0){
#      cat("overdispersion parameter =", fround (.Call("mer_sigma", 
#        object, FALSE, PACKAGE = "lme4"), 1), "\n")
#    }
#    }
#)


setMethod("display", signature(object = "mer"),
    function(object, digits=2, detail=FALSE)
    {
    call <- object@call
    print (call)
    #object <- summary(object)
    #summ <- summary(object)
    fcoef <- fixef(object)
    #coefs <- attr(summ, "coefs")
    #useScale <- attr (VarCorr (object), "sc")
    useScale <- object@dims["useSc"]
    corF <- vcov(object)@factors$correlation
    coefs <- cbind(fcoef, corF@sd)
    if (length (fcoef) > 0){
      if (!object@dims["useSc"]) {
        coefs <- coefs[, 1:2, drop = FALSE]
        stat <- coefs[, 1]/coefs[, 2]
        pval <- 2 * pnorm(abs(stat), lower = FALSE)
        coefs <- cbind(coefs, `z value` = stat, `Pr(>|z|)` = pval)
      }
      else {
        stat <- coefs[, 1]/coefs[, 2]
        coefs <- cbind(coefs, `t value` = stat)
      }
    dimnames(coefs)[[2]][1:2] <- c("coef.est", "coef.se")
      if(detail){
        pfround (coefs, digits)
      }
      else{
        pfround(coefs[,1:2], digits)
      }
    }
    cat("\nError terms:\n")
    vc <- as.matrix.VarCorr (VarCorr (object), useScale=useScale, digits)
    print (vc[,c(1:2,4:ncol(vc))], quote=FALSE)
    ngrps <- lapply(object@flist, function(x) length(levels(x)))
    REML <- object@dims["REML"]
    llik <- logLik(object, REML)
    AIC <- AIC(llik)
    dev <- object@deviance["ML"]     # Dbar
    n <- object@dims["n"]
    Dhat <- -2*(llik) # Dhat
    pD <- dev - Dhat              # pD
    DIC <- dev + pD               # DIC=Dbar+pD=Dhat+2pD
    cat("---\n")
    cat(sprintf("number of obs: %d, groups: ", n))
    cat(paste(paste(names(ngrps), ngrps, sep = ", "), collapse = "; "))
    cat(sprintf("\nAIC = %g, DIC = ", round(AIC,1)))
    cat(round(DIC, 1))
    cat("\ndeviance =", fround (dev, 1), "\n")
    if (useScale < 0){
      cat("overdispersion parameter =", fround (.Call("mer_sigma", 
        object, FALSE, PACKAGE = "lme4"), 1), "\n")
    }
    }
)


setMethod("display", signature(object = "polr"),
    function(object, digits=2, detail=FALSE)
    {
    call <- object$call
    summ <- summary(object)
    if(detail){
      coef <- summ$coef[, , drop = FALSE]
    }
    else{
      coef <- summ$coef[, 1:2, drop = FALSE]
    }
    dimnames(coef)[[2]][1:2] <- c("coef.est", "coef.se")
    n <- summ$n  
    k <- nrow (coef)
    k.intercepts <- length (summ$zeta)
    print(call)
    pfround(coef, digits)
    cat("---\n")
    cat(paste("n = ", n, ", k = ", k, " (including ", k.intercepts,
        " intercepts)\nresidual deviance = ",
        fround(deviance(object), 1), 
        ", null deviance is not computed by polr",
        "\n", sep = ""))
    #cat("AIC:", fround(AIC(object), 1), "\n")
    }
)


#setMethod("display", signature(object = "bayespolr"),
#    function(object, digits=2)
#    {
#    call <- object$call
#    summ <- summary(object)
#    coef <- summ$coef[, 1:2, drop = FALSE]
#    dimnames(coef)[[2]] <- c("coef.est", "coef.se")
#    n <- summ$n  # or maybe should be "nobs", I don't know for sure
#    k <- nrow (coef)
#    k.intercepts <- length (summ$zeta)
#    print(call)
#    pfround(coef, digits)
#    cat("---\n")
#    cat(paste("n = ", n, ", k = ", k, " (including ", k.intercepts,
#        " intercepts)\nresidual deviance = ",
#        fround(summ$deviance, 1), 
#        ", null deviance is not computed by bayespolr",
#        "\n", sep = ""))
#    }
#)
