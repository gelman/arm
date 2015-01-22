#' Extract AIC and DIC from a mer model
#'
#' Computes the (generalized) Akaike *A*n *I*nformation *C*riterion
#' and *D*eviance *I*nformation *C*riterion for a mer model.
#'
#' @param fit fitted \code{merMod} mode, usually the result of a
#'   fitter like \code{merMod}.
#' @param ... further arguments (currently unused).
#'
#' @author
#'   Andrew Gelman \email{gelman@@stat.columbia.edu};
#'   Yu-Sung Su \email{suyusung@@tsinghua.edu.cn}
#' @keywords manip methods
#' @export
#' @examples
#' fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
#' extractAIC(fm1)
#' extractDIC(fm1)
extractDIC <- function(fit,...){
  UseMethod("extractDIC")
}

#' @rdname extractDIC
#' @export
extractDIC.merMod <- function(fit,...){
#        REML <- fit@dims["REML"]
#        llik <- logLik(fit, REML)
#        dev <- fit@deviance["ML"]
#        n <- fit@dims["n"]
#        Dhat <- -2 * (llik)
#        pD <- dev - Dhat
#        DIC <- dev + pD[[1]]
#        names(DIC) <- "DIC"
#        return(DIC)      
        is_REML <- isREML(fit)
        llik <- logLik(fit, REML=is_REML)
        dev <- deviance(refitML(fit))
        n <-  getME(fit, "devcomp")$dims["n"]
        Dhat <- -2 * (llik)
        pD <- dev - Dhat
        DIC <- dev + pD[[1]]
        names(DIC) <- "DIC"
        return(DIC)
}




# #' @rdname extractDIC
# #' @export
#extractAIC.mer <- function(fit,...){
##     REML <- fit@dims["REML"]
##    llik <- logLik(fit, REML)
##    AIC <- AIC(llik)
##    names(AIC) <- "AIC"
##    return(AIC)
#    L <- logLik(refitML(fit))
#    edf <- attr(L,"df")
#    out <- c(edf,-2*L + k*edf)
#    return(out)
#} 
