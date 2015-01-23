## New bayespolr() using Kenny's Dirichlet prior distribution

#' Bayesian Ordered Logistic or Probit Regression
#'
#' Bayesian functions for ordered logistic or probit modeling with
#' independent normal, t, or Cauchy prior distribution for the
#' coefficients.
#'
#' The program is a simple alteration of \code{\link[MASS]{polr}} in
#' \code{VR} version 7.2-31 that augments the loglikelihood with the
#' log of the t prior distributions for the coefficients.
#'
#' We use Student-t prior distributions for the coefficients. The
#' prior distributions for the intercepts (the cutpoints) are set so
#' they apply to the value when all predictors are set to their mean
#' values.
#'
#' If scaled=TRUE, the scales for the prior distributions of the
#' coefficients are determined as follows: For a predictor with only
#' one value, we just use \code{prior.scale}. For a predictor with two
#' values, we use prior.scale/range(x). For a predictor with more than
#' two values, we use prior.scale/(2*sd(x)).
#'
#' @param formula a formula expression as for regression models, of
#'   the form \code{response ~ predictors}. The response should be a
#'   factor (preferably an ordered factor), which will be interpreted as
#'   an ordinal response, with levels ordered as in the factor.  A
#'   proportional odds model will be fitted. The model must have an
#'   intercept: attempts to remove one will lead to a warning and be
#'   ignored. An offset may be used.  See the documentation of
#'   \code{formula} for other details.
#' @param data an optional data frame in which to interpret the
#'   variables occurring in \code{formula}.
#' @param weights optional case weights in fitting.  Default to 1.
#' @param start initial values for the parameters.  This is in the
#'   format \code{c(coefficients, zeta)}
#' @param \dots additional arguments to be passed to \code{optim},
#'   most often a \code{control} argument.
#' @param subset expression saying which subset of the rows of the
#'   data should be used in the fit.  All observations are included by
#'   default.
#' @param na.action a function to filter missing data.
#' @param contrasts a list of contrasts to be used for some or all of
#'   the factors appearing as variables in the model formula.
#' @param Hess logical for whether the Hessian (the observed
#'   information matrix) should be returned.
#' @param model logical for whether the model matrix should be
#'   returned.
#' @param method logistic or probit or complementary log-log or
#'   cauchit (corresponding to a Cauchy latent variable and only
#'   available in R >= 2.1.0).
#' @param drop.unused.levels default \code{TRUE}, if \code{FALSE}, it
#'   interpolates the intermediate values if the data have integer
#'   levels.
#' @param prior.mean prior mean for the coefficients: default is 0.
#'   Can be a vector of length equal to the number of predictors (not
#'   counting the intercepts). If it is a scalar, it is expanded to the
#'   length of this vector.
#' @param prior.scale prior scale for the coefficients: default is
#'   2.5. Can be a vector of length equal to the number of predictors
#'   (not counting the intercepts). If it is a scalar, it is expanded to
#'   the length of this vector.
#' @param prior.df for t distribution: default is 1 (Cauchy).  Set to
#'   \code{Inf} to get normal prior distributions. Can be a vector of
#'   length equal to the number of predictors (not counting the
#'   intercepts). If it is a scalar, it is expanded to the length of
#'   this vector.
#' @param prior.counts.for.bins default is \code{NULL}, which will
#'   augment the data by giving each cut point a \code{1/levels(y)}. To
#'   use a noninformative prior, assign prior.counts.for.bins = 0. If it
#'   is a scalar, it is expanded to the number of levels of y.
#' @param min.prior.scale Minimum prior scale for the coefficients:
#'   default is 1e-12.
#' @param scaled if \code{scaled = TRUE}, then the prior distribution
#'   is rescaled.  Can be a vector of length equal to the number of
#'   cutpoints (intercepts). If it is a scalar, it is expanded to the
#'   length of this vector.
#' @param n.iter default is 100.
#' @param print.unnormalized.log.posterior display the unnormalized
#'   log posterior likelihood for bayesglm fit, default=\code{FALSE}
#'
#' @return See \code{polr} for details.
#'   \item{prior.mean}{prior means for the cofficients.}
#'   \item{prior.scale}{prior scales for the cofficients.}
#'   \item{prior.df}{prior dfs for the cofficients.}
#'   \item{prior.counts.for.bins}{prior counts for the cutpoints.}
#' @author
#'   Andrew Gelman \email{gelman@@stat.columbia.edu};
#'   Yu-Sung Su \email{suyusung@@tsinghua.edu.cn};
#'   Maria Grazia Pittau \email{grazia@@stat.columbia.edu}
#' @seealso \code{\link{bayesglm}}, \code{\link[MASS]{polr}}
#' @keywords models methods regression
#' @export
#' @examples
#' library("MASS")
#' M1 <- polr(Sat ~ Infl + Type + Cont, weights = Freq, data = housing)
#' display (M1)
#' 
#' M2 <- bayespolr(Sat ~ Infl + Type + Cont, weights = Freq, data = housing,
#'     prior.scale=Inf, prior.df=Inf) # Same as M1
#' display (M2)
#' 
#' M3 <- bayespolr(Sat ~ Infl + Type + Cont, weights = Freq, data = housing)
#' display (M3)
#' 
#' M4 <- bayespolr(Sat ~ Infl + Type + Cont, weights = Freq, data = housing,
#'     prior.scale=2.5, prior.df=1)  # Same as M3
#' display (M4)   
#' 
#' M5 <- bayespolr(Sat ~ Infl + Type + Cont, weights = Freq, data = housing,
#'     prior.scale=2.5, prior.df=7)
#' display (M5)
#'
#' M6 <- bayespolr(Sat ~ Infl + Type + Cont, weights = Freq, data = housing,
#'     prior.scale=2.5, prior.df=Inf)
#' display (M6)
#' 
#' # Assign priors 
#' M7 <- bayespolr(Sat ~ Infl + Type + Cont, weights = Freq, data = housing,
#'     prior.mean=rep(0,6), prior.scale=rep(2.5,6), prior.df=c(1,1,1,7,7,7))
#' display (M7)
#' 
#' 
#' #### Another example
#' y <- factor (rep (1:10,1:10))
#' x <- rnorm (length(y))
#' x <- x - mean(x)
#' 
#' M8 <- polr (y ~ x)
#' display (M8)
#'
#' M9 <- bayespolr (y ~ x,  prior.scale=Inf, prior.df=Inf, prior.counts.for.bins=0)   
#' display (M9) # same as M1
#'
#' M10 <- bayespolr (y ~ x,  prior.scale=Inf, prior.df=Inf, prior.counts.for.bins=10000)
#' display (M10)
#'
#'
#' #### Another example
#'
#' y <- factor (rep (1:3,1:3))
#' x <- rnorm (length(y))
#' x <- x - mean(x)
#'
#' M11 <- polr (y ~ x)
#' display (M11)
#'
#' M12 <- bayespolr (y ~ x,  prior.scale=Inf, prior.df=Inf, prior.counts.for.bins=0)   
#' display (M12) # same as M1
#'
#' M13 <- bayespolr (y ~ x,  prior.scale=Inf, prior.df=Inf, prior.counts.for.bins=1)
#' display (M13)
#'
#' M14 <- bayespolr (y ~ x,  prior.scale=Inf, prior.df=Inf, prior.counts.for.bins=10)
#' display (M14)
bayespolr <- 
function (formula, data, weights, start, ..., subset, na.action, 
    contrasts = NULL, Hess = TRUE, model = TRUE, method = c("logistic", 
        "probit", "cloglog", "cauchit"), drop.unused.levels = TRUE, 
    prior.mean = 0, prior.scale = 2.5, prior.df = 1, prior.counts.for.bins = NULL,
    min.prior.scale = 1e-12,
    scaled = TRUE, n.iter = 100, print.unnormalized.log.posterior = FALSE) 
{
    logit <- function(p) log(p/(1 - p))

    dt.deriv <- function(x, mean, scale, df, log = TRUE, delta = 0.001) {
        (dt((x + delta - mean)/scale, df, log = log) - dt((x - 
            delta - mean)/scale, df, log = log))/(2 * delta)
    }

    fmin <- function(beta) {
        theta <- beta[pc + 1:q]
        gamm <- c(-100, cumsum(c(theta[1], exp(theta[-1]))), 
            100)
        eta <- offset
        if (pc > 0)
            eta <- eta + drop(x %*% beta[1:pc])
        pr <- pfun(gamm[y + 1] - eta) - pfun(gamm[y] - eta)
        if (all(pr > 0)) 
            f <- -sum(wt * log(pr))
        else f <- Inf
        if (pc > 0) 
            f <- f - sum(dt((beta[1:pc] - prior.mean)/prior.scale, 
                prior.df, log = TRUE))
        return(f)
    }

    gmin <- function(beta) {
        jacobian <- function(theta) {
            k <- length(theta)
            etheta <- exp(theta)
            mat <- matrix(0, k, k)
            mat[, 1] <- rep(1, k)
            for (i in 2:k) mat[i:k, i] <- etheta[i]
            mat
        }
        theta <- beta[pc + 1:q]
        gamm <- c(-100, cumsum(c(theta[1], exp(theta[-1]))), 
            100)
        eta <- offset
        if (pc > 0) 
            eta <- eta + drop(x %*% beta[1:pc])
        pr <- pfun(gamm[y + 1] - eta) - pfun(gamm[y] - eta)
        p1 <- dfun(gamm[y + 1] - eta)
        p2 <- dfun(gamm[y] - eta)
        g1 <- if (pc > 0) 
            t(x) %*% (wt * (p1 - p2)/pr)
        else numeric(0)
        xx <- .polrY1 * p1 - .polrY2 * p2
        g2 <- -t(xx) %*% (wt/pr)
        g2 <- t(g2) %*% jacobian(theta)
        if (pc > 0)
            g1 <- g1 - dt.deriv(beta[1:pc], prior.mean, prior.scale, 
                prior.df, log = TRUE)
        if (all(pr > 0)) 
            c(g1, g2)
        else rep(NA, pc + q)
    }
    m <- match.call(expand.dots = FALSE)
    mf <- match(c("formula", "data", "subset", "weights", "na.action", 
        "etastart", "mustart", "offset"), names(m), 0)
    m <- m[c(1, mf)]
    m$drop.unused.levels <- drop.unused.levels
    method <- match.arg(method)
  
  ##### adjust prior.scale for probit ####
  if (method == "probit"){
    prior.scale <- prior.scale*1.6
  }
  ################
  
  for(jj in 1:length(prior.scale)){
    if (prior.scale[jj] < min.prior.scale){
      prior.scale[jj] <- min.prior.scale
      warning ("prior scale for variable ", jj, " set to min.prior.scale = ", min.prior.scale,"\n")
    }
  }

  
  

    pfun <- switch(method, logistic = plogis, probit = pnorm, 
        cloglog = pgumbel, cauchit = pcauchy)
    dfun <- switch(method, logistic = dlogis, probit = dnorm, 
        cloglog = dgumbel, cauchit = dcauchy)
    if (is.matrix(eval.parent(m$data))) 
        m$data <- as.data.frame(data)
    m$start <- m$Hess <- m$method <- m$... <- NULL
    m[[1]] <- as.name("model.frame")
    m <- eval.parent(m)
    Terms <- attr(m, "terms")
    x <- model.matrix(Terms, m, contrasts)
    xint <- match("(Intercept)", colnames(x), nomatch = 0)
    n <- nrow(x)
    pc <- ncol(x)
    cons <- attr(x, "contrasts")
    if (xint > 0) {
        x <- x[, -xint, drop = FALSE]
        pc <- pc - 1
    }
    else warning("an intercept is needed and assumed")
    wt <- model.weights(m)
    if (!length(wt)) 
        wt <- rep(1, n)
    offset <- model.offset(m)
    if (length(offset) <= 1) 
        offset <- rep(0, n)
    y <- model.response(m)
    if (!is.factor(y)) 
        stop("response must be a factor")
    lev <- levels(y)
    if (length(lev) <= 2) 
        stop("response must have 3 or more levels")
    y <- unclass(y)
    q <- length(lev) - 1
    Y <- matrix(0, n, q)
    .polrY1 <- col(Y) == y
    .polrY2 <- col(Y) == y - 1
    if (missing(start)) {
        q1 <- length(lev)%/%2
        y1 <- (y > q1)
        X <- cbind(Intercept = rep(1, n), x)
        fit <- switch(method, 
            logistic = bayesglm.fit(X, y1, 
              wt, family = binomial(), offset = offset, intercept = TRUE, 
              prior.mean = prior.mean, prior.scale = prior.scale, 
              prior.df = prior.df, prior.mean.for.intercept = 0, 
              prior.scale.for.intercept = 10, prior.df.for.intercept = 1, 
              min.prior.scale = min.prior.scale,
              scaled = scaled, control = glm.control(maxit = n.iter), 
              print.unnormalized.log.posterior = print.unnormalized.log.posterior), 
            probit = bayesglm.fit(X, y1, wt, family = binomial("probit"), 
                offset = offset, intercept = TRUE, prior.mean = prior.mean, 
                prior.scale = prior.scale, prior.df = prior.df, 
                prior.mean.for.intercept = 0, prior.scale.for.intercept = 10, 
                prior.df.for.intercept = 1, 
                min.prior.scale = min.prior.scale,
                scaled = scaled, 
                control = glm.control(maxit = n.iter), print.unnormalized.log.posterior = print.unnormalized.log.posterior), 
            cloglog = bayesglm.fit(X, y1, wt, family = binomial("probit"), 
                offset = offset, intercept = TRUE, prior.mean = prior.mean, 
                prior.scale = prior.scale, prior.df = prior.df, 
                prior.mean.for.intercept = 0, prior.scale.for.intercept = 10, 
                prior.df.for.intercept = 1, 
                min.prior.scale = min.prior.scale,
                scaled = scaled, 
                control = glm.control(maxit = n.iter), print.unnormalized.log.posterior = print.unnormalized.log.posterior), 
            cauchit = bayesglm.fit(X, y1, wt, family = binomial("cauchit"), 
                offset = offset, intercept = TRUE, prior.mean = prior.mean, 
                prior.scale = prior.scale, prior.df = prior.df, 
                prior.mean.for.intercept = 0, prior.scale.for.intercept = 10, 
                prior.df.for.intercept = 1, 
                min.prior.scale = min.prior.scale,
                scaled = scaled, 
                control = glm.control(maxit = n.iter), print.unnormalized.log.posterior = print.unnormalized.log.posterior))
        if (!fit$converged) 
            warning("attempt to find suitable starting values failed")
        coefs <- fit$coefficients
        if (any(is.na(coefs))) {
            warning("design appears to be rank-deficient, so dropping some coefs")
            keep <- names(coefs)[!is.na(coefs)]
            coefs <- coefs[keep]
            x <- x[, keep[-1], drop = FALSE]
            pc <- ncol(x)
        }
        spacing <- logit((1:q)/(q + 1))
        if (method != "logistic") 
            spacing <- spacing/1.7
        gammas <- -coefs[1] + spacing - spacing[q1]
        thetas <- c(gammas[1], log(diff(gammas)))
        start <- c(coefs[-1], thetas)
    }
    # rep start to have the same length of coef + zeta
    else if (length(start)==1){
      start <- rep(start, (pc+q))
    }
    else if (length(start) != pc + q) 
        stop("'start' is not of the correct length")
    
    J <- NCOL(x)
    
    # SU: if no x's, no priors for coefs  2008.2.9
    if (xint>1) {
      if (length(prior.mean) == 1) 
        prior.mean <- rep(prior.mean, J)
      if (length(prior.scale) == 1) {
        prior.scale <- rep(prior.scale, J)
          if (scaled == TRUE) {
            for (j in 1:J) {
              x.obs <- x[, j]
              x.obs <- x.obs[!is.na(x.obs)]
              num.categories <- length(unique(x.obs))
              if (num.categories == 2) {
                prior.scale[j] <- prior.scale[j]/(max(x.obs) - min(x.obs))
              }
              else if (num.categories > 2) {
                prior.scale[j] <- prior.scale[j]/(2 * sd(x.obs))
              }
            }
          }
      }
      if (length(prior.df) == 1) {
        prior.df <- rep(prior.df, J)
      }
    }

    # prior for intercept sum(priors.intercpet)=1
    if (is.null(prior.counts.for.bins)) {
      prior.counts.for.bins <- 1/(q+1)
    }
    if (length(prior.counts.for.bins) == 1) {
      prior.counts.for.bins <- rep(prior.counts.for.bins, q+1)
    }
# Augment the data to add prior information
    y.0 <- y
    Y.0 <- Y
    x.0 <- x
    wt.0 <- wt
    offset.0 <- offset
    .polrY1.0 <- .polrY1
    .polrY2.0 <- .polrY2
    y <- c (y.0, 1:(q+1))
    Y <- matrix(0, n+q+1, q)
    .polrY1 <- col(Y) == y
    .polrY2 <- col(Y) == y - 1
    x <- rbind (x.0, matrix (colMeans(x.0), nrow=(q+1), ncol=J, byrow=TRUE))
    wt <- c (wt.0, prior.counts.for.bins)
    offset <- c (offset, rep(0,q+1))
# Fit the model as before
    res <- optim(start, fmin, gmin, method = "BFGS", hessian = Hess, ...)
# Restore the old variables
    y <- y.0
    Y <- Y.0
    x <- x.0
    wt <- wt.0
    offset <- offset.0
    .polrY1 <- .polrY1.0
    .polrY2 <- .polrY2.0
# Continue on as before
    beta <- res$par[seq_len(pc)]
    theta <- res$par[pc + 1:q]
    zeta <- cumsum(c(theta[1], exp(theta[-1])))
    deviance <- 2 * res$value
    niter <- c(f.evals = res$counts[1], g.evals = res$counts[2])
    names(zeta) <- paste(lev[-length(lev)], lev[-1], sep = "|")
    if (pc > 0) {
        names(beta) <- colnames(x)
        eta <- drop(x %*% beta)
    }
    else {
        eta <- rep(0, n)
    }
    cumpr <- matrix(pfun(matrix(zeta, n, q, byrow = TRUE) - eta), , q)
    fitted <- t(apply(cumpr, 1, function(x) diff(c(0, x, 1))))
    dimnames(fitted) <- list(row.names(m), lev)
    fit <- list(coefficients = beta, zeta = zeta, deviance = deviance, 
        fitted.values = fitted, lev = lev, terms = Terms, df.residual = sum(wt) - 
            pc - q, edf = pc + q, n = sum(wt), nobs = sum(wt), 
        call = match.call(), method = method, convergence = res$convergence, 
        prior.mean = prior.mean, prior.scale = prior.scale, prior.df = prior.df, 
        prior.counts.for.bins = prior.counts.for.bins, niter = niter)
    if (Hess) {
        dn <- c(names(beta), names(zeta))
        H <- res$hessian
        dimnames(H) <- list(dn, dn)
        fit$Hessian <- H
    }
    if (model){
        fit$model <- m
    }
    fit$na.action <- attr(m, "na.action")
    fit$contrasts <- cons
    fit$xlevels <- .getXlevels(Terms, m)
    class(fit) <- c("bayespolr", "polr")
    fit
}

#' @rdname bayespolr
#' @export
setMethod("print", signature(x = "bayespolr"), 
    function(x, digits= 2) display(object=x, digits=digits))

#' @rdname bayespolr
#' @export
setMethod("show", signature(object = "bayespolr"), 
    function(object) display(object, digits=2))
