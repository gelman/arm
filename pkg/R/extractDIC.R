extractDIC <- function(fit,...){
  UseMethod("extractDIC")
}


extractDIC.mer <- function(fit,...){
        REML <- fit@dims["REML"]
        llik <- logLik(fit, REML)
        dev <- fit@deviance["ML"]
        n <- fit@dims["n"]
        Dhat <- -2 * (llik)
        pD <- dev - Dhat
        DIC <- dev + pD[[1]]
        names(DIC) <- "DIC"
        return(DIC)
}





extractAIC.mer <- function(fit,...){
        REML <- fit@dims["REML"]
        llik <- logLik(fit, REML)
        AIC <- AIC(llik)
        names(AIC) <- "AIC"
        return(AIC)
} 
