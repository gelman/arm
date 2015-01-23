
#' Functions for Processing lm, glm, mer, polr and svyglm Output
#'
#' This generic function gives a clean printout of lm, glm, mer, polr
#' and svyglm objects.
#'
#' This generic function gives a clean printout of lm, glm, mer and
#' polr objects, focusing on the most pertinent pieces of information:
#' the coefficients and their standard errors, the sample size, number
#' of predictors, residual standard deviation, and R-squared. Note:
#' R-squared is automatically displayed to 2 digits, and deviances are
#' automatically displayed to 1 digit, no matter what.
#'
#' @note Output are the model, the regression coefficients and
#' standard errors, and the residual sd and R-squared (for a linear
#' model), or the null deviance and residual deviance (for a
#' generalized linear model).
#'
#' @param object The output of a call to lm, glm, mer, polr, svyglm or
#'   related regressions function with n data points and k predictors.
#' @param ... further arguments passed to or from other methods.
#' @param digits number of significant digits to display.
#' @param detail defaul is \code{FALSE}, if \code{TRUE}, display
#'   p-values or z-values.
#' 
#' @return Coefficients and their standard errors, the sample size,
#'   number of predictors, residual standard deviation, and R-squared
#' @references Andrew Gelman and Jennifer Hill, Data Analysis Using
#'   Regression and Multilevel/Hierarchical Models, Cambridge University
#'   Press, 2006.
#' @author
#'   Andrew Gelman \email{gelman@@stat.columbia.edu}; 
#'   Yu-Sung Su \email{suyusung@@tsinghua.edu.cn};
#'   Maria Grazia Pittau \email{grazia@@stat.columbia.edu}
#' @keywords manip methods
#' @name display
#' @examples
#' # Here's a simple example of a model of the form, y = a + bx + error, 
#' # with 10 observations in each of 10 groups, and with both the 
#' # intercept and the slope varying by group.  First we set up the model and data.
#'   library("MASS")  # Needed for mvrnorm() and polr()
#'   group <- rep(1:10, rep(10,10))
#'   group2 <- rep(1:10, 10)
#'   mu.a <- 0
#'   sigma.a <- 2
#'   mu.b <- 3
#'   sigma.b <- 4
#'   rho <- 0.56
#'   Sigma.ab <- array (c(sigma.a^2, rho*sigma.a*sigma.b, 
#'                    rho*sigma.a*sigma.b, sigma.b^2), c(2,2))
#'   sigma.y <- 1
#'   ab <- mvrnorm (10, c(mu.a,mu.b), Sigma.ab)
#'   a <- ab[,1]
#'   b <- ab[,2]
#'   d <- rnorm(10)
#'
#'   x <- rnorm (100)
#'   y1 <- rnorm (100, a[group] + b*x, sigma.y)
#'   y2 <- rbinom(100, 1, prob=invlogit(a[group] + b*x))
#'   y3 <- rnorm (100, a[group] + b[group]*x + d[group2], sigma.y)
#'   y4 <- rbinom(100, 1, prob=invlogit(a[group] + b*x + d[group2]))
#'
#'
#' # display a simple linear model
#' 
#'    M1 <- lm (y1 ~ x)
#'    display (M1)
#'    M1.sim <- sim(M1, n.sims=2)  
#' 
#' # display a simple logit model
#' 
#'    M2 <- glm (y2 ~ x, family=binomial(link="logit"))
#'    display (M2)
#'    M2.sim <- sim(M2, n.sims=2)  
#' 
#' # Then fit and display a simple varying-intercept model:
#'  
#'    library("lme4")
#'    M3 <- lmer (y1 ~ x + (1|group))
#'    display (M3)
#'    M3.sim <- sim(M3, n.sims=2)  
#' 
#' 
#' # Then the full varying-intercept, varying-slope model:
#'  
#'    M4 <- lmer (y1 ~ x + (1 + x |group))
#'    display (M4)
#'    M4.sim <- sim(M4, n.sims=2)  
#' 
#'  
#' # Then the full varying-intercept, logit model:
#'  
#'    M5 <- glmer (y2 ~ x + (1|group), family=binomial(link="logit"))
#'    display (M5)
#'    M5.sim <- sim(M5, n.sims=2)  
#' 
#'  
#' # Then the full varying-intercept, varying-slope logit model:
#'  
#'    M6 <- glmer (y2 ~ x + (1|group) + (0 + x |group), 
#'         family=binomial(link="logit"))
#'    display (M6)
#'    M6.sim <- sim(M6, n.sims=2)  
#' 
#' 
#' # Then non-nested varying-intercept, varying-slop model:
#'    
#'    M7 <- lmer (y3 ~ x + (1 + x |group) + (1|group2))
#'    display(M7)
#'    M7.sim <- sim(M7, n.sims=2)  
#' 
#'  
#' # Then the ordered logit model from polr
#' 
#'    M8 <- polr(Sat ~ Infl + Type + Cont, weights = Freq, data = housing)
#'    display(M8)
#'    
#'    M9 <- bayespolr(Sat ~ Infl + Type + Cont, weights = Freq, data = housing)
#'    display(M9) 
NULL


#' @rdname display
#' @export
setMethod("display", signature(object = "lm"),
    function(object, digits=2, detail=FALSE)
    {
    out <- NULL
    out$call <- object$call
    summ <- summary (object)
    out$sigma.hat <- summ$sigma
    out$r.squared <- summ$r.squared
    if(detail){
      coef <- summ$coef[,,drop=FALSE]
    }
    else{
      coef <- summ$coef[,1:2,drop=FALSE]
    }
    dimnames(coef)[[2]][1:2] <- c("coef.est","coef.se")
    out$coef <- coef[,"coef.est"]#,drop=FALSE]
    out$se <- coef[,"coef.se"]#,drop=FALSE]
    out$t.value <- summ$coef[,3]
    out$p.value <- summ$coef[,4]
    out$n <- summ$df[1] + summ$df[2]
    out$k <- summ$df[1]
    print (out$call)
    pfround (coef, digits)
    cat("---\n")
    cat (paste ("n = ", out$n, ", k = ", out$k,
    "\nresidual sd = ", fround (out$sigma.hat, digits),
    ", R-Squared = ", fround (out$r.squared, 2), "\n", sep=""))
    return(invisible(out))
    }
)


#' @rdname display
#' @export
setMethod("display", signature(object = "bayesglm"),
    function(object, digits=2, detail=FALSE)
    {
    out <- NULL
    out$call <- object$call
    summ <- summary(object, dispersion = object$dispersion)
    if(detail){
      coef <- summ$coefficients
      coef[ rownames( coef ) %in% rownames( summ$coef[, , drop = FALSE]) , ] <- summ$coef[ , , drop = FALSE ] 
      out$z.value <- coef[,3]#,drop=FALSE]
      out$p.value <- coef[,4]#,drop=FALSE]
    }
    else{
      coef <- matrix( NA, length( object$coefficients ),2 )
      rownames(coef) <- names( object$coefficients )          ## M
      coef[ rownames( coef ) %in% rownames( summ$coef[, 1:2, drop = FALSE]) , ] <- summ$coef[ , 1:2, drop = FALSE ]  ## M
    }
    dimnames(coef)[[2]][1:2] <- c( "coef.est", "coef.se")
    out$coef <- coef[,"coef.est"]#,drop=FALSE]
    out$se <- coef[,"coef.se"]#,drop=FALSE]
    out$n <- summ$df[1] + summ$df[2]
    out$k <- summ$df[1]
    out$deviance <- summ$deviance
    out$null.deviance <- summ$null.deviance
    print(out$call)
    pfround(coef, digits)
    cat("---\n")
    cat(paste("n = ", out$n, ", k = ", out$k, "\nresidual deviance = ", 
        fround(out$deviance, 1), ", null deviance = ", fround(out$null.deviance, 1), " (difference = ", fround(out$null.deviance - out$deviance, 1), ")", "\n", sep = ""))
    out$dispersion <- if (is.null(object$dispersion)){
                        summ$dispersion
                      } else {
                        object$dispersion
                      }
    if (out$dispersion != 1) {
        out$overdispersion.parameter <- out$dispersion
        cat(paste("overdispersion parameter = ", fround(out$dispersion, 1), "\n", sep = ""))
        if (family(object)$family == "gaussian") {
          out$sigma.hat <- sqrt(out$dispersion)
          cat(paste("residual sd is sqrt(overdispersion) = ", fround(out$sigma.hat, digits), "\n", sep = ""))
        }
    }
    return(invisible(out))
    }
)


#setMethod("display", signature(object = "bayesglm.h"),
#    function (object, digits = 2, detail = FALSE) 
#    {
#    call <- object$call
#    summ <- summary(object, dispersion = object$dispersion)
#    if(detail){
#      coef <- summ$coefficients
#      coef[ rownames( coef ) %in% rownames( summ$coef[, , drop = FALSE]) , ] <- summ$coef[ , , drop = FALSE ] 
#    }
#    else{
#      coef <- matrix( NA, length( object$coefficients ),2 )
#      rownames(coef) <- names( object$coefficients )          ## M
#      coef[ rownames( coef ) %in% rownames( summ$coef[, 1:2, drop = FALSE]) , ] <- summ$coef[ , 1:2, drop = FALSE ]  ## M
#    }
#    dimnames(coef)[[2]][1:2] <- c( "coef.est", "coef.se")
#    #n <- summ$df[1] + summ$df[2]
#    n <- summ$df.residual
#    k <- summ$df[1]
#    print(call)
#    if(max(object$batch)>0){
#        nn<- strsplit( rownames( coef )[seq( from= length( object$batch ) + 1 ,to = nrow( coef ))], "." , fixed=TRUE)
#        bb<- c( object$batch,unlist( lapply (nn , function( lst ) { lst[[3]] } ) ) ) 
#    }
#    else {bb<- c( object$batch)}
#    cc<- cbind( fround( coef, digits ), bb )
#    dimnames(cc)[[2]][3]<-"batch"
#    print( cc , quote = FALSE )
#    cat("---\n")
#    cat(paste("n = ", n, ", k = ", k, "\nresidual deviance = ", 
#        fround(summ$deviance, 1), ", null deviance = ", fround(summ$null.deviance, 
#            1), " (difference = ", fround(summ$null.deviance - 
#            summ$deviance, 1), ")", "\n", sep = ""))
#    dispersion <- if (is.null(object$dispersion)) 
#        summ$dispersion
#    else object$dispersion
#    if (dispersion != 1) {
#        cat(paste("overdispersion parameter = ", fround(dispersion, 
#            1), "\n", sep = ""))
#        if (family(object)$family == "gaussian") {
#            cat(paste("residual sd is sqrt(overdispersion) = ", 
#                fround(sqrt(dispersion), digits), "\n", sep = ""))
#            cat(paste("group sd is sigma.batch = ", 
#                fround(object$sigma.batch, digits), "\n", sep = ""))
#        }
#    }
#    }
#)


#' @rdname display
#' @export
setMethod("display", signature(object = "glm"),
    function(object, digits=2, detail=FALSE)
    {
    out <- NULL
    out$call <- object$call
    summ <- summary(object, dispersion = object$dispersion)
    if(detail){
      coef <- summ$coef[, , drop = FALSE]
      out$z.value <- coef[,3]#,drop=FALSE]
      out$p.value <- coef[,4]#,drop=FALSE]
    }
    else{
      coef <- summ$coef[, 1:2, drop = FALSE]
    }
    dimnames(coef)[[2]][1:2] <- c("coef.est", "coef.se")
    out$n <- summ$df[1] + summ$df[2]
    out$k <- summ$df[1]
    out$coef <- coef[,"coef.est"]
    out$se <- coef[,"coef.se"]
     print(out$call)
    pfround(coef, digits)
    out$deviance <- summ$deviance
    out$null.deviance <- summ$null.deviance
    cat("---\n")
    cat(paste("  n = ", out$n, ", k = ", out$k, "\n  residual deviance = ", 
        fround(out$deviance, 1), ", null deviance = ", fround(out$null.deviance, 1), " (difference = ", fround(summ$null.deviance - summ$deviance, 1), ")", "\n", sep = ""))
    out$dispersion <- if (is.null(object$dispersion)){
                        summ$dispersion
                      } else {
                        object$dispersion
                      }
    if (out$dispersion != 1) {
      cat(paste("  overdispersion parameter = ", fround(out$dispersion, 1), "\n", sep = ""))
      if (family(object)$family=="gaussian") {
        out$sigma.hat <- sqrt(out$dispersion)
        cat(paste("  residual sd is sqrt(overdispersion) = ",
                  fround(out$sigma.hat, digits), "\n", sep = ""))
      }
    }
    return(invisible(out))
  }
)




#setMethod("display", signature(object = "mer"),
#    function(object, digits=2)
#    {
#    call <- object@call
#    print (call)
#    #object <- summary(object)
#    fcoef <- fixef(object)
#    useScale <- attr( lme4::VarCorr(object), "sc") 
#    corF <- vcov(object)@factors$correlation
#    coefs <- cbind(fcoef, corF@sd)
#    if (length (fcoef) > 0){
#      dimnames(coefs) <- list(names(fcoef), c("coef.est", "coef.se"))
#      pfround (coefs, digits)
#    }
#    cat("\nError terms:\n")
#    vc <- as.matrix.VarCorr (lme4::VarCorr (object), useScale=useScale, digits)
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



#' @rdname display
#' @export
setMethod("display", signature(object = "merMod"),
    function(object, digits=2, detail=FALSE)
    {
    out <- NULL
    out$call <- object@call
    print (out$call)
    #object <- summary(object)
    #summ <- summary(object)
    fcoef <- fixef(object)
    #coefs <- attr(summ, "coefs")
    #useScale <- attr (lme4::VarCorr (object), "sc")
    useScale <- getME(object, "devcomp")$dims["useSc"]
    corF <- vcov(object)@factors$correlation
    coefs <- cbind(fcoef, corF@sd)
    if (length (fcoef) > 0){
      if (!useScale) {
        coefs <- coefs[, 1:2, drop = FALSE]
        out$z.value <- coefs[, 1]/coefs[, 2]
        out$p.value <- 2 * pnorm(abs(out$z.value), lower = FALSE)
        coefs <- cbind(coefs, `z value` = out$z.value, `Pr(>|z|)` = out$p.value)
      }
      else {
        out$t.value <- coefs[, 1]/coefs[, 2]
        coefs <- cbind(coefs, `t value` = out$t.value)
      }
    dimnames(coefs)[[2]][1:2] <- c("coef.est", "coef.se")
      if(detail){
        pfround (coefs, digits)
      }
      else{
        pfround(coefs[,1:2], digits)
      }
    }
    out$coef <- coefs[,"coef.est"]
    out$se <- coefs[,"coef.se"]
    cat("\nError terms:\n")
    vc <- as.matrix.VarCorr (lme4::VarCorr (object), useScale=useScale, digits)
    print (vc[,c(1:2,4:ncol(vc))], quote=FALSE)
    out$ngrps <- lapply(object@flist, function(x) length(levels(x)))
    is_REML <- isREML(object)
    llik <- logLik(object, REML=is_REML)
    out$AIC <- AIC(llik)
    out$deviance <- deviance(refitML(object))     # Dbar
    out$n <- getME(object, "devcomp")$dims["n"]
    Dhat <- -2*(llik) # Dhat
    pD <- out$deviance - Dhat              # pD
    out$DIC <- out$deviance + pD               # DIC=Dbar+pD=Dhat+2pD
    cat("---\n")
    cat(sprintf("number of obs: %d, groups: ", out$n))
    cat(paste(paste(names(out$ngrps), out$ngrps, sep = ", "), collapse = "; "))
    cat(sprintf("\nAIC = %g, DIC = ", round(out$AIC,1)))
    cat(round(out$DIC, 1))
    cat("\ndeviance =", fround (out$deviance, 1), "\n")
    if (useScale < 0){
      out$sigma.hat <- .Call("mer_sigma", object, FALSE, PACKAGE = "lme4")
      cat("overdispersion parameter =", fround (out$sigma.hat, 1), "\n")
    }
    return(invisible(out))
    }
)


#' @rdname display
#' @export
setMethod("display", signature(object = "polr"),
    function(object, digits=2, detail=FALSE)
    {
    out <- NULL
    out$call <- object$call
    summ <- summary(object)
    if(detail){
      coef <- summ$coef[, , drop = FALSE]
      out$t.value <- coef[,"t value"]
    }
    else{
      coef <- summ$coef[, 1:2, drop = FALSE]
    }
    dimnames(coef)[[2]][1:2] <- c("coef.est", "coef.se")
    out$coef <- coef[,"coef.est"]
    out$se <- coef[,"coef.se"]
    out$n <- summ$n  
    out$k <- nrow (coef)
    out$k.intercepts <- length (summ$zeta)
    print(out$call)
    pfround(coef, digits)
    cat("---\n")
    cat(paste("n = ", out$n, ", k = ", out$k, " (including ", out$k.intercepts,
        " intercepts)\nresidual deviance = ",
        fround(deviance(object), 1), 
        ", null deviance is not computed by polr",
        "\n", sep = ""))
    #cat("AIC:", fround(AIC(object), 1), "\n")
    return(invisible(out))
    }
)


#' @rdname display
#' @export
setMethod("display", signature(object = "svyglm"),
    function(object, digits=2, detail=FALSE)
    {
    out <- NULL
    out$call <- object$call
    out$survey.design <- object$survey.design
    summ <- summary(object)
    if(detail){
      coef <- summ$coef[, , drop = FALSE]
      out$z.value <- coef[,3]#,drop=FALSE]
      out$p.value <- coef[,4]#,drop=FALSE]
    }
    else{
      coef <- summ$coef[, 1:2, drop = FALSE]
    }
    dimnames(coef)[[2]][1:2] <- c("coef.est", "coef.se")
    out$n <- summ$df[1] + summ$df[2]
    out$k <- summ$df[1]
    out$coef <- coef[,"coef.est"]
    out$se <- coef[,"coef.se"]
    print(out$call)
    cat("\n")
    print(out$survey.design)
    cat("\n")
    pfround(coef, digits)
    out$deviance <- summ$deviance
    out$null.deviance <- summ$null.deviance
    cat("---\n")
    cat(paste("  n = ", out$n, ", k = ", out$k, "\n  residual deviance = ", 
        fround(out$deviance, 1), ", null deviance = ", fround(out$null.deviance, 1), " (difference = ", fround(summ$null.deviance - summ$deviance, 1), ")", "\n", sep = ""))
    out$dispersion <- summ$dispersion[1]
    if (out$dispersion != 1) {
      cat(paste("  overdispersion parameter = ", fround(out$dispersion, 1), "\n", sep = ""))
      if (family(object)$family=="gaussian") {
        out$sigma.hat <- sqrt(out$dispersion)
        cat(paste("  residual sd is sqrt(overdispersion) = ",
                  fround(out$sigma.hat, digits), "\n", sep = ""))
      }
    }
    return(invisible(out))
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
