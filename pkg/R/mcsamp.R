#' Generic Function to Run mcmcsamp() in lme4
#'
#' Quick function for MCMC sampling for lmer and glmer objects and
#' convert to Bugs objects for easy display.
#'
#' This function generates a sample from the posterior distribution of
#' the parameters of a fitted model using Markov Chain Monte Carlo
#' methods. It automatically simulates multiple sequences and allows
#' convergence to be monitored. The function relies on
#' \code{\link[lme4]{mcmcsamp}} in \code{lme4}.
#'
#' @param object \code{mer} objects from \code{lme4}
#' @param n.chains number of MCMC chains.
#' @param n.iter number of iteration for each MCMC chain.
#' @param n.burnin number of burnin for each MCMC chain, Default is
#' \code{n.iter/2}, that is, discarding the first half of the
#' simulations.
#' @param n.thin keep every kth draw from each MCMC chain. Must be a
#' positive integer. Default is \code{max(1, floor(n.chains *
#' (n.iter-n.burnin) / 1000))} which will only thin if there are at
#' least 2000 simulations.
#' @param saveb if 'TRUE', causes the values of the random effects in
#' each sample to be saved.
#' @param deviance compute deviance for \code{mer} objects. Only works
#' for \code{\link[lme4]{lmer}} object.
#' @param make.bugs.object tranform the output into bugs object,
#' default is TRUE.
#' @param ... further arguments passed to or from other methods.
#' @return An object of (S3) class '"bugs"' suitable for use with the
#' functions in the "R2WinBUGS" package.
#' @references Andrew Gelman and Jennifer Hill, Data Analysis Using
#' Regression and Multilevel/Hierarchical Models, Cambridge University
#' Press, 2006.
#'
#' Douglas Bates and Deepayan Sarkar, lme4: Linear mixed-effects
#' models using S4 classes.
#' @author
#'   Andrew Gelman \email{gelman@@stat.columbia.edu};
#'   Yu-Sung Su \email{ys463@@columbia.edu}
#' @seealso \code{\link{display}}, \code{\link[lme4]{lmer}},
#' \code{\link[lme4]{mcmcsamp}}, \code{\link{sim}}
#' @keywords models methods
#' @name mcsamp
#' @export
#' @examples
#' \dontrun{
#' # Here's a simple example of a model of the form, y = a + bx + error, 
#' # with 10 observations in each of 10 groups, and with both the intercept 
#' # and the slope varying by group.  First we set up the model and data.
#' #   
#' group <- rep(1:10, rep(10,10))
#' group2 <- rep(1:10, 10)
#' mu.a <- 0
#' sigma.a <- 2
#' mu.b <- 3
#' sigma.b <- 4
#' rho <- 0.56
#' Sigma.ab <- array (c(sigma.a^2, rho*sigma.a*sigma.b, 
#'                  rho*sigma.a*sigma.b, sigma.b^2), c(2,2))
#' sigma.y <- 1
#' ab <- mvrnorm (10, c(mu.a,mu.b), Sigma.ab)
#' a <- ab[,1]
#' b <- ab[,2]
#' d <- rnorm(10)
#' 
#' x <- rnorm (100)
#' y1 <- rnorm (100, a[group] + b*x, sigma.y)
#' y2 <- rbinom(100, 1, prob=invlogit(a[group] + b*x))
#' y3 <- rnorm (100, a[group] + b[group]*x + d[group2], sigma.y)
#' y4 <- rbinom(100, 1, prob=invlogit(a[group] + b*x + d[group2]))
#' 
#' # 
#' # Then fit and display a simple varying-intercept model:
#'  
#' M1 <- lmer (y1 ~ x + (1|group))
#' display (M1)
#' M1.sim <- mcsamp (M1)
#' print (M1.sim)
#' plot (M1.sim)
#' 
#' # 
#' # Then the full varying-intercept, varying-slope model:
#' # 
#' M2 <- lmer (y1 ~ x + (1 + x |group))
#' display (M2)
#' M2.sim <- mcsamp (M2)
#' print (M2.sim)
#' plot (M2.sim)
#' 
#' # 
#' # Then the full varying-intercept, logit model:
#' # 
#' M3 <- lmer (y2 ~ x + (1|group), family=binomial(link="logit"))
#' display (M3)
#' M3.sim <- mcsamp (M3)
#' print (M3.sim)
#' plot (M3.sim)
#' 
#' # 
#' # Then the full varying-intercept, varying-slope logit model:
#' # 
#' M4 <- lmer (y2 ~ x + (1|group) + (0+x |group), 
#'      family=binomial(link="logit"))
#' display (M4)
#' M4.sim <- mcsamp (M4)
#' print (M4.sim)
#' plot (M4.sim)
#'    
#' #
#' # Then non-nested varying-intercept, varying-slop model:
#' #
#' M5 <- lmer (y3 ~ x + (1 + x |group) + (1|group2))
#' display(M5)
#' M5.sim <- mcsamp (M5)
#' print (M5.sim)
#' plot (M5.sim)
#' }
NULL

mcsamp.default <- function (object, n.chains=3, n.iter=1000, n.burnin=floor(n.iter/2), 
    n.thin=max(1, floor(n.chains * (n.iter - n.burnin)/1000)), 
    saveb=TRUE, deviance=TRUE, make.bugs.object=TRUE)
{
  cat("mcsamp() used to be a wrapper for mcmcsamp() in lme4.\nCurrently, mcmcsamp() is no longer available in lme4.\nSo in the meantime, we suggest that users use sim() to get\nsimulated estimates.\n")
}



#mcsamp.default <- function (object, n.chains=3, n.iter=1000, n.burnin=floor(n.iter/2), 
#    n.thin=max(1, floor(n.chains * (n.iter - n.burnin)/1000)), 
#    saveb=TRUE, deviance=TRUE, make.bugs.object=TRUE)
#{
#  
#  if (n.chains<2) stop ("n.chains must be at least 2")
#  n.keep <- n.iter - n.burnin
#  first.chain <- mcmcsamp (object, n.iter, saveb=saveb, trans=TRUE, deviance=deviance)[(n.burnin+1):n.iter,]
#  n.parameters <- ncol(first.chain)
#    
#  if (deviance) {
#    sims <- array (NA, c(n.keep, n.chains, n.parameters+1))
#  }
#  if (!deviance){
#    sims <- array (NA, c(n.keep, n.chains, n.parameters))
#  }
#
#  pred.names <- attr(terms(object), "term.labels")
#  par.names <- dimnames(first.chain)[[2]]
#  par.names <- gsub("b.", "b@", par.names, ignore.case = FALSE, # Su: rename "b.*" to ""
#                    extended = TRUE, perl = FALSE,
#                    fixed = TRUE, useBytes = FALSE)   
#  par.names <- gsub("b@.*", "", par.names, ignore.case = FALSE, 
#                    extended = TRUE, perl = FALSE,
#                    fixed = FALSE)    
#  par.names <- par.names[is.na(match(par.names,""))] 
#  name.chk.idx <- as.logical(match(par.names, pred.names, nomatch=0))
#  par.names[name.chk.idx] <- paste("beta", par.names[name.chk.idx], sep=".")
#
#  if (saveb){
#    b.hat <- se.coef (object)                   # Su: use se.coef() 
#    n.groupings <- length(b.hat) - 1
#    J <- NA
#    K <- NA
#    for (m in 1:n.groupings){
#      J[m] <- dim(b.hat[[m+1]])[1]
#      K[m] <- dim(b.hat[[m+1]])[2]
#      var.names <- paste (abbreviate(names(b.hat)[m+1],4), ".",
#                          unlist (dimnames(b.hat[[m+1]])[2]), sep="") ##sep="."
#      par.names <- c (par.names,
#        paste (rep(var.names,J[m]), "[", rep(1:J[m],each=K[m]), "]", sep=""))
#    }
#  }
#  sims[,1,1:n.parameters] <- first.chain
#
#  for (k in 2:n.chains){
#    sims[,k,1:n.parameters] <- mcmcsamp (object, n.iter, saveb=saveb, trans=TRUE, deviance=deviance)[(n.burnin+1):n.iter,]
#  }
#  
#  select <- c(rep(FALSE, n.thin-1),TRUE)
#  sims <- sims[select,,]
#  
#  for (j in 1:n.parameters){
#    if (pmatch("log(sigma^2)", par.names[j], nomatch=0)){#=="log(sigma^2)"){
#      par.names[j] <- "sigma.y"
#      sims[,,j] <- exp (sims[,,j]/2)
#    }
#    else if (pmatch("log(", par.names[j], nomatch=0)){#(substr(par.names[j],1,4)=="log("){
#      par.names[j] <- paste ("sigma.", substr(par.names[j], 5, nchar(par.names[j])-1), sep="")
#      sims[,,j] <- exp (sims[,,j]/2)
#    }
#    else if (pmatch("atanh(", par.names[j], nomatch=0)){#(substr(par.names[j],1,6)=="atanh("){
#      par.names[j] <- paste ("rho.", substr(par.names[j], 7, nchar(par.names[j])-1), sep="")
#      sims[,,j] <- tanh (sims[,,j])
#    }
#    #else if (substr(par.names[j],1,4)=="eta."){#(pmatch("eta.", par.names[j], nomatch=0)){#(substr(par.names[j],1,4)=="eta."){
#    #  par.names[j] <- paste ("", substr(par.names[j], 5, nchar(par.names[j])), sep="")
#    #  par.names[j] <- par.names[j]
#    #}
#    else if (pmatch("deviance", par.names[j], nomatch=0)){#(par.names[j]=="deviance"){          # Su: keep par.names for "deviance"
#        sims[,,n.parameters+1] <- sims[,,j] 
#        sims <- sims[,,-j]                          # Su: delete deviance value from sims
#    } 
##    else {  
##    }
#  } 
#  par.names <- gsub("(", "", par.names, ignore.case = FALSE,
#                    extended = TRUE, perl = FALSE,
#                    fixed = TRUE, useBytes = FALSE)   
#  par.names <- gsub(")", "", par.names, ignore.case = FALSE,
#                    extended = TRUE, perl = FALSE,
#                    fixed = TRUE, useBytes = FALSE) 
# # par.names <- gsub(".Intercept", ".Int", par.names, ignore.case = FALSE,
##                    extended = TRUE, perl = FALSE,
##                    fixed = TRUE, useBytes = FALSE) 
#  par.names <- gsub("rescale", "z.", par.names, ignore.case = FALSE,
#                    extended = TRUE, perl = FALSE,
#                    fixed = TRUE, useBytes = FALSE) 
#  
#  par.names <- par.names[is.na(match(par.names,"deviance"))] # Su: delete par.names for "deviance"
#  
#  if (deviance){
#      dimnames(sims) <- list (NULL, NULL, c(par.names,"deviance"))
#  }
#  if (!deviance){
#    dimnames(sims) <- list (NULL, NULL, par.names)
#  }
#  if (make.bugs.object){
#    return (as.bugs.array (sims, program="lmer", n.iter=n.iter, n.burnin=n.burnin, n.thin=n.thin, DIC=deviance))
#  }
#  else {
#    return (sims)
#  }
#}
#
#
#

#' @rdname mcsamp
#' @export
setMethod("mcsamp", signature(object = "merMod"),
    function (object, n.chains=3, n.iter=1000, n.burnin=floor(n.iter/2), 
              n.thin=max(1, floor(n.chains * (n.iter - n.burnin)/1000)), 
              saveb=TRUE, deviance=TRUE, make.bugs.object=TRUE)
{
    mcsamp.default(object, deviance=TRUE, ...)
}
)
#
#setMethod("mcsamp", signature(object = "glmer"),
#    function (object, ...)
#{
#    mcsamp.default(object, deviance=FALSE, ...)
#}
#)
