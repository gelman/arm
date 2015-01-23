#' Functions to Get Posterior Distributions
#'
#' This generic function gets posterior simulations of sigma and beta from a
#' \code{lm} object, or simulations of beta from a \code{glm} object, or
#' simulations of beta from a \code{merMod} object
#'
#' @param object the output of a call to \code{lm} with n data points
#'   and k predictors.
#' @param slot return which slot of \code{sim.polr}, available options
#'   are \code{coef, zeta, ALL}.
#' @param ... further arguments passed to or from other methods.
#' @param n.sims number of independent simulation draws to create.
#' @param regression the orginial mer model.
#' @return
#'   \item{coef}{matrix (dimensions n.sims x k) of n.sims random
#'     draws of coefficients.}
#'   \item{zeta}{matrix (dimensions n.sims x k) of n.sims random draws of
#'     zetas (cut points in polr).}
#'   \item{fixef}{matrix (dimensions n.sims x k) of n.sims random draws
#'     of coefficients of the fixed effects for the \code{merMod}
#'     objects. Previously, it is called \code{unmodeled}.}
#'   \item{sigma}{vector of n.sims random draws of sigma (for \code{glm}'s, this
#'     just returns a vector of 1's or else of the square root of the
#'     overdispersion parameter if that is in the model)}
#' @references Andrew Gelman and Jennifer Hill. (2006). \emph{Data
#'   Analysis Using Regression and Multilevel/Hierarchical
#'   Models}. Cambridge University Press.
#' @author
#'   Andrew Gelman \email{gelman@@stat.columbia.edu};
#'   Yu-Sung Su \email{suyusung@@tsinghua.edu.cn};
#'   Vincent Dorie \email{vjd4@@nyu.edu}
#' @seealso \code{\link{display}}, \code{\link{lm}},
#' \code{\link{glm}}, \code{\link[lme4]{lmer}}
#' @keywords models methods
#' @name sim
#' @examples
#' #Examples of "sim" 
#'  set.seed (1)
#'  J <- 15
#'  n <- J*(J+1)/2
#'  group <- rep (1:J, 1:J)
#'  mu.a <- 5
#'  sigma.a <- 2
#'  a <- rnorm (J, mu.a, sigma.a)
#'  b <- -3
#'  x <- rnorm (n, 2, 1)
#'  sigma.y <- 6
#'  y <- rnorm (n, a[group] + b*x, sigma.y)
#'  u <- runif (J, 0, 3)
#'  y123.dat <- cbind (y, x, group)
#' 
#' # Linear regression 
#'  x1 <- y123.dat[,2]
#'  y1 <- y123.dat[,1]
#'  M1 <- lm (y1 ~ x1)
#'  display(M1)
#'  M1.sim <- sim(M1)
#'  coef.M1.sim <- coef(M1.sim)
#'  sigma.M1.sim <- sigma.hat(M1.sim)
#'  ## to get the uncertainty for the simulated estimates
#'  apply(coef(M1.sim), 2, quantile)
#'  quantile(sigma.hat(M1.sim))
#'  
#' # Logistic regression 
#'  u.data <- cbind (1:J, u)
#'  dimnames(u.data)[[2]] <- c("group", "u")
#'  u.dat <- as.data.frame (u.data)
#'  y <- rbinom (n, 1, invlogit (a[group] + b*x))
#'  M2 <- glm (y ~ x, family=binomial(link="logit"))
#'  display(M2)
#'  M2.sim <- sim (M2)
#'  coef.M2.sim <- coef(M2.sim)
#'  sigma.M2.sim <- sigma.hat(M2.sim)
#' 
#' # Ordered Logistic regression 
#'  library("MASS")
#'  house.plr <- polr(Sat ~ Infl + Type + Cont, weights = Freq, data = housing)
#'  display(house.plr)
#'  M.plr <- sim(house.plr)
#'  coef.sim <- coef(M.plr, slot="coef")   
#'  zeta.sim <- coef(M.plr, slot="zeta")
#'  coefall.sim <- coef(M.plr)
#' 
#' # Using lmer:
#' # Example 1
#'  library("lme4")
#'  E1 <- lmer (y ~ x + (1 | group))
#'  display(E1)
#'  E1.sim <- sim (E1)
#'  coef.E1.sim <- coef(E1.sim)
#'  fixef.E1.sim <- fixef(E1.sim)
#'  ranef.E1.sim <- ranef(E1.sim)
#'  sigma.E1.sim <- sigma.hat(E1.sim)
#'  yhat <- fitted(E1.sim, E1)
#' 
#' # Example 2
#'  u.full <- u[group]
#'  E2 <- lmer (y ~ x + u.full + (1 | group))
#'  display(E2)
#'  E2.sim <- sim (E2)
#'  coef.E2.sim <- coef(E2.sim)
#'  fixef.E2.sim <- fixef(E2.sim)
#'  ranef.E2.sim <- ranef(E2.sim)
#'  sigma.E2.sim <- sigma.hat(E2.sim)
#'  yhat <- fitted(E2.sim, E2)
#' 
#' # Example 3 
#'  y <- rbinom (n, 1, invlogit (a[group] + b*x))
#'  E3 <- glmer (y ~ x + (1 | group), family=binomial(link="logit"))
#'  display(E3)
#'  E3.sim <- sim (E3)
#'  coef.E3.sim <- coef(E3.sim)
#'  fixef.E3.sim <- fixef(E3.sim)
#'  ranef.E3.sim <- ranef(E3.sim)
#'  sigma.E3.sim <- sigma.hat(E3.sim)
#'  yhat <- fitted(E3.sim, E3)
NULL

#' @rdname sim
#' @export
setMethod("sim", signature(object = "lm"),
    function(object, n.sims=100)
    {
    object.class <- class(object)[[1]]
    summ <- summary (object)
    coef <- summ$coef[,1:2,drop=FALSE]
    dimnames(coef)[[2]] <- c("coef.est","coef.sd")
    sigma.hat <- summ$sigma
    beta.hat <- coef[,1,drop = FALSE]
    V.beta <- summ$cov.unscaled
    n <- summ$df[1] + summ$df[2]
    k <- summ$df[1]
    sigma <- rep (NA, n.sims)
    beta <- array (NA, c(n.sims,k))
    dimnames(beta) <- list (NULL, rownames(beta.hat))
    for (s in 1:n.sims){
      sigma[s] <- sigma.hat*sqrt((n-k)/rchisq(1,n-k))
      beta[s,] <- mvrnorm (1, beta.hat, V.beta*sigma[s]^2)
    }

    ans <- new("sim",
                coef = beta,
                sigma = sigma)
    return (ans)
    }
)

#' @rdname sim
#' @export
setMethod("sim", signature(object = "glm"),
    function(object, n.sims=100)
    {
    object.class <- class(object)[[1]]
    summ <- summary (object, correlation=TRUE, dispersion = object$dispersion)
    coef <- summ$coef[,1:2,drop=FALSE]
    dimnames(coef)[[2]] <- c("coef.est","coef.sd")
    beta.hat <- coef[,1,drop=FALSE]
    sd.beta <- coef[,2,drop=FALSE]
    corr.beta <- summ$corr
    n <- summ$df[1] + summ$df[2]
    k <- summ$df[1]
    V.beta <- corr.beta * array(sd.beta,c(k,k)) * t(array(sd.beta,c(k,k)))
    beta <- array (NA, c(n.sims,k))
    dimnames(beta) <- list (NULL, dimnames(beta.hat)[[1]])
    for (s in 1:n.sims){
      beta[s,] <- mvrnorm (1, beta.hat, V.beta)
    }
    # Added by Masanao
    beta2 <- array (0, c(n.sims,length(coefficients(object))))
    dimnames(beta2) <- list (NULL, names(coefficients(object)))
    beta2[,dimnames(beta2)[[2]]%in%dimnames(beta)[[2]]] <- beta
    # Added by Masanao
    sigma <- rep (sqrt(summ$dispersion), n.sims)

    ans <- new("sim",
                coef = beta2,
                sigma = sigma)
    return(ans)
    }
)

#' @rdname sim
#' @export
setMethod("sim", signature(object = "polr"),
    function(object, n.sims=100){
  x <- as.matrix(model.matrix(object))
  coefs <- coef(object)
  k <- length(coefs)
  zeta <- object$zeta
  Sigma <- vcov(object)

  if(n.sims==1){
    parameters <- t(mvrnorm(n.sims, c(coefs, zeta), Sigma))
  }else{
    parameters <- mvrnorm(n.sims, c(coefs, zeta), Sigma)
  }
  ans <- new("sim.polr",
              coef = parameters[,1:k,drop=FALSE],
              zeta = parameters[,-(1:k),drop=FALSE])
  return(ans)
})


#setMethod("sim", signature(object = "mer"),
#    function(object, n.sims=100)
#    {
#    #object <- summary(object)
##    if (lapply(object@bVar,sum)<=0|sum(unlist(lapply(object@bVar, is.na)))>0){
##        object@call$control <- list(usePQL=TRUE)
##        object <- lmer(object@call$formula)
#    #}
#    #sc <- attr (VarCorr (object), "sc")
#    # simulate unmodeled coefficients
#
#    fcoef <- fixef(object)
#    corF <- vcov(object)@factors$correlation
#    se.unmodeled <- corF@sd
#    V.beta <- (se.unmodeled %o% se.unmodeled) * as.matrix(corF)
#    beta.unmodeled <- NULL
#    if (length (fcoef) > 0){
#      beta.unmodeled[[1]] <- mvrnorm (n.sims, fcoef, V.beta)
#      names (beta.unmodeled) <- "unmodeled"
#    }
#    # simulate coefficients within groups
#    #coef <- ranef (object)
#    #estimate <- ranef(object, postVar=TRUE)
#    #vars <- object@bVar
#    #beta.bygroup <- vars
#
#    sc <- attr (VarCorr (object), "sc")
#    coef <- ranef(object, postVar=TRUE)
#    beta.bygroup <- c(coef)
#    n.groupings <- length (coef)
#    for (m in 1:n.groupings){
#      #vars.m <- vars[[m]]
#      vars.m <- attr (coef[[m]], "postVar")
#      K <- dim(vars.m)[1]
#      J <- dim(vars.m)[3]
#      beta.bygroup[[m]] <- array (NA, c(n.sims, J, K))
#      bhat <- coef[[m]]
#      for (j in 1:J){
#        V.beta <- untriangle(vars.m[,,j])#*sc^2
#        beta.bygroup[[m]][,j,] <- mvrnorm (n.sims, bhat[j,], V.beta)
#      }
#      dimnames (beta.bygroup[[m]]) <- c (list(NULL), dimnames(bhat))
#    }
#    betas <- c (beta.unmodeled, beta.bygroup)
#    return (betas)
#    }
#)

#setMethod("sim", signature(object = "mer"),
#    function(object, n.sims=100, ranef=TRUE)
#    {
#    # simulate unmodeled coefficients
#    fcoef <- fixef(object)
#    corF <- vcov(object)@factors$correlation
#    se.unmodeled <- corF@sd
#    V.beta <- (se.unmodeled %o% se.unmodeled) * as.matrix(corF)
#    beta.unmodeled <- NULL
#    if (length (fcoef) > 0){
#      beta.unmodeled[[1]] <- mvrnorm (n.sims, fcoef, V.beta)
#      names (beta.unmodeled) <- "fixef"#"unmodeled"
#      coef <- beta.unmodeled
#    }
#    if(ranef){
#      # simulate coefficients within groups
#      sc <- attr (VarCorr (object), "sc")  # scale
#      #coef <- ranef (object)
#      #estimate <- ranef(object, postVar=TRUE)
#      coef <- ranef(object, postVar=TRUE)
#      beta.bygroup <- coef
#      n.groupings <- length (coef)
#      for (m in 1:n.groupings){
#        bhat <- as.matrix(coef[[m]]) # to suit the use of mvrnorm
#        vars.m <- attr (coef[[m]], "postVar")
#        K <- dim(vars.m)[1]
#        J <- dim(vars.m)[3]
#        beta.bygroup[[m]] <- array (NA, c(n.sims, J, K))
#        for (j in 1:J){
#          V.beta <- .untriangle(vars.m[,,j])#*sc^2
#          beta.bygroup[[m]][,j,] <- mvrnorm (n.sims, bhat[j,], V.beta)
#        }
#        dimnames (beta.bygroup[[m]]) <- c (list(NULL), dimnames(bhat))
#      }
#      coef <- c (beta.unmodeled, beta.bygroup)
#      }
#    return (coef)
#    }
#)
