bayesglm <- function (formula, family = gaussian, data, weights, subset,
  na.action, start = NULL, etastart, mustart, offset, control =
  glm.control(...),
  model = TRUE, method = "glm.fit", x = FALSE, y = TRUE, contrasts = NULL,
  drop.unused.levels = TRUE,
  prior.mean = 0, prior.scale = NULL, prior.df = 1,
  prior.mean.for.intercept = 0,
  prior.scale.for.intercept = NULL, prior.df.for.intercept = 1,
  min.prior.scale = 1e-12, scaled = TRUE, keep.order = TRUE,
  drop.baseline = TRUE, n.iter = 100, 
  print.unnormalized.log.posterior = FALSE, Warning=TRUE,...)
{
  call <- match.call()
  if (is.character(family))
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family))
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
  if (missing(data))
    data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action", 
    "etastart", "mustart", "offset"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- drop.unused.levels
  mf$na.action <- NULL
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  switch(method, model.frame = return(mf), glm.fit = 1, stop("invalid 'method': ",
    method))
  mt <- attr(mf, "terms")
  Y <- model.response(mf, "any")
  if (length(dim(Y)) == 1) {
    nm <- rownames(Y)
    dim(Y) <- NULL
    if (!is.null(nm))
      names(Y) <- nm
  }
  X <- if (!is.empty.model(mt)) {
#      class(mt) <- c("bayesglm", oldClass(mt))
      model.matrix.bayes(mt, mf, contrasts, keep.order = keep.order, drop.baseline=drop.baseline)
    }
    else matrix(, NROW(Y), 0)
#  }
#  else {
#    X <- if (!is.empty.model(mt))
#      model.matrix(mt, mf, contrasts)
#  else matrix(, NROW(Y), 0)
#  }
  weights <- model.weights(mf)
  offset <- model.offset(mf)
  if (!is.null(weights) && any(weights < 0))
    stop("negative weights not allowed")
  if (!is.null(offset) && length(offset) != NROW(Y))
    stop(gettextf("number of offsets is %d should equal %d (number of observations)",
      length(offset), NROW(Y)), domain = NA)
  mustart <- model.extract(mf, "mustart")
  etastart <- model.extract(mf, "etastart")
  fit <- bayesglm.fit(x = X, y = Y, weights = weights, start = start,
    etastart = etastart, mustart = mustart, offset = offset,
    family = family, control = glm.control(maxit = n.iter),
    intercept = attr(mt, "intercept") > 0, prior.mean = prior.mean,
    prior.scale = prior.scale, prior.df = prior.df, 
    prior.mean.for.intercept = prior.mean.for.intercept,
    prior.scale.for.intercept = prior.scale.for.intercept,
    prior.df.for.intercept = prior.df.for.intercept, min.prior.scale =
    min.prior.scale, scaled = scaled, Warning=Warning)
  if (any(offset) && attr(mt, "intercept") > 0) {
   cat("bayesglm not yet set up to do deviance comparion here\n")
  fit$null.deviance <- bayesglm.fit(x = X[, "(Intercept)", drop = FALSE], 
    y = Y, weights = weights, offset = offset,
    family = family, control = control, intercept = TRUE,
    prior.mean = prior.mean, prior.scale = prior.scale,
    prior.df = prior.df, prior.mean.for.intercept = prior.mean.for.intercept,
    prior.scale.for.intercept = prior.scale.for.intercept,
    prior.df.for.intercept = prior.df.for.intercept,
    min.prior.scale = min.prior.scale, scaled = scaled, 
    Warning=Warning)$deviance
  }
  if (model)
    fit$model <- mf
  fit$na.action <- attr(mf, "na.action")
  if (x)
    fit$x <- X
  if (!y)
    fit$y <- NULL
  fit <- c(fit, list(call = call, formula = formula, terms = mt,
    data = data, offset = offset, control = control, method = method,
    contrasts = attr(X, "contrasts"), xlevels = .getXlevels(mt, mf)))
  class(fit) <- c("bayesglm", "glm", "lm")
  fit
  }
  
bayesglm.fit <- function (x, y, weights = rep(1, nobs), start = NULL, 
  etastart = NULL, mustart = NULL, offset = rep(0, nobs), 
  family = gaussian(),
  control = glm.control(), intercept = TRUE,
  prior.mean = 0,
  prior.scale = NULL,
  prior.df = 1,
  prior.mean.for.intercept = 0,
  prior.scale.for.intercept = NULL,
  prior.df.for.intercept = 1,
  min.prior.scale=1e-12,
  scaled = TRUE, print.unnormalized.log.posterior=FALSE, Warning=TRUE)
{

  ##### 12.13 ####
  if (is.null(prior.scale)){
    prior.scale <- 2.5
  if (family$link == "probit")
    prior.scale <- prior.scale*1.6
  }
  #prior.scale <- prior.scale
  
  if (is.null(prior.scale.for.intercept)){
    prior.scale.for.intercept <- 10
  if (family$link == "probit")
    prior.scale.for.intercept <- prior.scale.for.intercept*1.6
  }
  #prior.scale.for.intercept <- prior.scale.for.intercept
  ################
  
  J <- NCOL(x)
  if (length(prior.mean) == 1) {
    prior.mean <- rep(prior.mean, J)
  }
  else {
    if (intercept) {
      prior.mean <- c(prior.mean.for.intercept, prior.mean)
    }
#   else {
#      prior.mean <- prior.mean
#    }
  }
  if (length(prior.scale)==1){
    prior.scale <- rep(prior.scale, J)
  }
  else {
    if (intercept) {
      prior.scale <- c(prior.scale.for.intercept, prior.scale)
    }
#    else {
#      prior.scale <- prior.scale
#    }
  }
  if (length(prior.df) == 1) {
    prior.df <- rep(prior.df, J)
  }
  else {
    if (intercept) {
      prior.df <- c(prior.df.for.intercept, prior.df)
    }
#    else {
#      prior.df <- prior.df
#    }
  }
  
  
  if (scaled) {
    if (family$family == "gaussian")
      prior.scale <- prior.scale * 2 * sd(y)
    prior.scale.0 <- prior.scale
  
    for (j in 1:J) {
      x.obs <- x[, j]
      x.obs <- x.obs[!is.na(x.obs)]
      num.categories <- length(unique(x.obs))
      x.scale <- 1
      if (num.categories == 2) {
        x.scale <- max(x.obs) - min(x.obs)
      }
      else if (num.categories > 2) {
        x.scale <- 2 * sd(x.obs)
      }
      prior.scale[j] <- prior.scale[j]/x.scale
#      if(is.na(prior.scale[j])){
#        prior.scale[j] <- min.prior.scale
#        warning ("prior scale for variable ", j,
#          " set to min.prior.scale = ", min.prior.scale,"\n")      
#      }
      if (prior.scale[j] < min.prior.scale){
        prior.scale[j] <- min.prior.scale
        warning ("prior scale for variable ", j,
          " set to min.prior.scale = ", min.prior.scale,"\n")
       }
    }
  }
  
  x <- as.matrix(x)
  xnames <- dimnames(x)[[2]]
  ynames <- if (is.matrix(y)) rownames(y) else names(y)
  conv <- FALSE
  nobs <- NROW(y)
  nvars <- ncol(x)
  EMPTY <- nvars == 0
  if (is.null(weights))
    weights <- rep.int(1, nobs)
  if (is.null(offset))
    offset <- rep.int(0, nobs)
  variance <- family$variance
  dev.resids <- family$dev.resids
  aic <- family$aic
  linkinv <- family$linkinv
  mu.eta <- family$mu.eta
  if (!is.function(variance) || !is.function(linkinv))
    stop("'family' argument seems not to be a valid family object")
  valideta <- family$valideta
  if (is.null(valideta))
    valideta <- function(eta) TRUE
    validmu <- family$validmu
  if (is.null(validmu))
    validmu <- function(mu) TRUE
  if (is.null(mustart))
    eval(family$initialize)
  else {
    mukeep <- mustart
    eval(family$initialize)
    mustart <- mukeep
  }
  if (EMPTY) {
    eta <- rep.int(0, nobs) + offset
    if (!valideta(eta))
      stop("invalid linear predictor values in empty model")
    mu <- linkinv(eta)
    if (!validmu(mu))
      stop("invalid fitted means in empty model")
    dev <- sum(dev.resids(y, mu, weights))
    w <- ((weights * mu.eta(eta)^2)/variance(mu))^0.5
    residuals <- (y - mu)/mu.eta(eta)
    good <- rep(TRUE, length(residuals))
    boundary <- conv <- TRUE
    coef <- numeric(0)
    iter <- 0
  }
  else {
    coefold <- NULL
    if (!is.null(etastart)) {
      eta <- etastart
    }
    else if (!is.null(start)) {
      if (length(start) != nvars)
        stop(gettextf("length of 'start' should equal %d and correspond to initial coefs for %s", 
          nvars, paste(deparse(xnames), collapse = ", ")), domain = NA)
      else {
        eta <- offset + as.vector(ifelse((NCOL(x) == 1), x *start, x %*% start))
        coefold <- start
      }
    }
    else {
      eta <- family$linkfun(mustart)
    }
    mu <- linkinv(eta)
    if (!(validmu(mu) && valideta(eta)))
      stop("cannot find valid starting values: please specify some")
    devold <- sum(dev.resids(y, mu, weights))
    boundary <- conv <- FALSE
    prior.sd <- prior.scale
    #========Andy 2008.7.8=============
    dispersion <- ifelse((family$family %in% c("poisson", "binomial")),  1, var(y)/10000)
    #==================================
    dispersionold <- dispersion
    for (iter in 1:control$maxit) {
      good <- weights > 0
      varmu <- variance(mu)[good]
      if (any(is.na(varmu)))
        stop("NAs in V(mu)")
      if (any(varmu == 0))
        stop("0s in V(mu)")
      mu.eta.val <- mu.eta(eta)
      if (any(is.na(mu.eta.val[good])))
        stop("NAs in d(mu)/d(eta)")
      good <- (weights > 0) & (mu.eta.val != 0)
      if (all(!good)) {
        conv <- FALSE
        warning("no observations informative at iteration ", iter)
        break
      }
      z <- (eta - offset)[good] + (y - mu)[good]/mu.eta.val[good]
      w <- sqrt((weights[good] * mu.eta.val[good]^2)/variance(mu)[good])
      ngoodobs <- as.integer(nobs - sum(!good))
      coefs.hat <- rep(0, ncol(x))
      z.star <- c(z, prior.mean)
      x.extra <- diag(NCOL(x))
      if (intercept & scaled) {
        x.extra[1,] <- colMeans(x)
      }
      x.star <- rbind (x, x.extra)
      w.star <- c(w, sqrt(dispersion)/prior.sd)
      good.star <- c(good, rep(TRUE, NCOL(x)))
      ngoodobs.star <- ngoodobs + NCOL(x)
      fit <- .Fortran("dqrls", qr = x.star[good.star, ] * w.star, 
        n = ngoodobs.star, p = nvars, y = w.star * z.star, 
        ny = as.integer(1), tol = min(1e-07, control$epsilon/1000), 
        coefficients = double(nvars), residuals = double(ngoodobs.star), 
        effects = double(ngoodobs.star), rank = integer(1), 
        pivot = 1:nvars, qraux = double(nvars),
        work = double(2 * nvars), PACKAGE = "base")
      if (any(!is.finite(fit$coefficients))) {
        conv <- FALSE
        warning("non-finite coefficients at iteration ", iter)
        break
      }
      coefs.hat <- fit$coefficients
      fit$qr <- as.matrix (fit$qr)
      V.coefs <- chol2inv(fit$qr[1:ncol(x.star), 1:ncol(x.star), drop = FALSE])
      #V.beta <- chol2inv (t(x.star) %*% diag(w.star^2) %*% x.star)
  
      if (family$family == "gaussian" & scaled){
        prior.scale <- prior.scale.0 #*sqrt(dispersion)
      }
      centered.coefs <- coefs.hat
      sampling.var <- diag(V.coefs)
      if (intercept & scaled){
        centered.coefs[1] <- sum(coefs.hat*colMeans(x))
        sampling.var[1] <- t(colMeans(x))%*%V.coefs%*%colMeans(x)
      }
  
  # Andy 2007.12.13
      prior.sd <- ifelse(prior.df == Inf, prior.scale,
      sqrt(((centered.coefs - prior.mean)^2 + sampling.var * dispersion + 
        prior.df * prior.scale^2)/(1 + prior.df)))  
  # wrong????
      start[fit$pivot] <- fit$coefficients
      eta <- drop(x %*% start)
      mu <- linkinv(eta <- eta + offset)
      dev <- sum(dev.resids(y, mu, weights))
  
      if (!(family$family %in% c("poisson", "binomial"))) {
        mse.resid <- mean((w * (z - x %*% coefs.hat))^2)
  #mse.uncertainty <- mean(diag(x %*% V.coefs %*% t(x))) * dispersion
        mse.uncertainty <- mean(rowSums(( x %*% V.coefs ) * x)) * dispersion #faster
        dispersion <- mse.resid + ifelse(mse.uncertainty<0, 0, mse.uncertainty)
      }
  
      if (control$trace)
        cat("Deviance =", dev, "Iterations -", iter, "\n")
          boundary <- FALSE
  
      if (!is.finite(dev)) {
        if (is.null(coefold))
          stop("no valid set of coefficients has been found: please supply starting values", 
            call. = FALSE)
        warning("step size truncated due to divergence", call. = FALSE)
        ii <- 1
        while (!is.finite(dev)) {
          if (ii > control$maxit)
            stop("inner loop 1; cannot correct step size")
          ii <- ii + 1
          start <- (start + coefold)/2
          eta <- drop(x %*% start)
          mu <- linkinv(eta <- eta + offset)
          dev <- sum(dev.resids(y, mu, weights))
        }
        boundary <- TRUE
        if (control$trace)
          cat("Step halved: new deviance =", dev, "\n")
      }
      if (!(valideta(eta) && validmu(mu))) {
        if (is.null(coefold))
          stop("no valid set of coefficients has been found: please supply starting values",
            call. = FALSE)
        warning("step size truncated: out of bounds", call. = FALSE)
        ii <- 1
        while (!(valideta(eta) && validmu(mu))) {
          if (ii > control$maxit)
            stop("inner loop 2; cannot correct step size")
          ii <- ii + 1
          start <- (start + coefold)/2
          eta <- drop(x %*% start)
          mu <- linkinv(eta <- eta + offset)
      }
      boundary <- TRUE
      dev <- sum(dev.resids(y, mu, weights))
      if (control$trace)
        cat("Step halved: new deviance =", dev, "\n")
      }
  
  ############################
      if ( family$family == "binomial" && print.unnormalized.log.posterior) {
        logprior <- if( intercept ) {
          sum( dt( coefs.hat[-1], prior.df ,prior.mean,log = TRUE ) )
            + dt( coefs.hat[1], prior.df.for.intercept, prior.mean.for.intercept,
            log = TRUE )
        }
        else {
          sum( dt( coefs.hat, prior.df ,prior.mean, log = TRUE ) )
        }
        xb <- invlogit( x %*% coefs.hat )
        loglikelihood <- sum( log( c( xb[ y == 1 ], 1 - xb[ y == 0 ] ) ) )
        cat( "log prior: ", logprior, ", log likelihood: ", loglikelihood, ",
          unnormalized log posterior: ", loglikelihood +logprior, "\n" ,sep="")
      }
  ####################
  
      if (iter > 1 & abs(dev - devold)/(0.1 + abs(dev)) <  control$epsilon & 
        abs(dispersion - dispersionold)/(0.1 + abs(dispersion)) < control$epsilon) {
        conv <- TRUE
        coef <- start
        break
      }
      else {
        devold <- dev
        dispersionold <- dispersion
        coef <- coefold <- start
      }
    }
    if(Warning){
      if (!conv)
        warning("algorithm did not converge")
      if (boundary)
        warning("algorithm stopped at boundary value")
      eps <- 10 * .Machine$double.eps
      if (family$family == "binomial") {
        if (any(mu > 1 - eps) || any(mu < eps)) {
          warning("fitted probabilities numerically 0 or 1 occurred")
        }
      }
      if (family$family == "poisson") {
        if (any(mu < eps))
          warning("fitted rates numerically 0 occurred")
      }
    }
    if (fit$rank < nvars) {
      coef[fit$pivot][seq(fit$rank + 1, nvars)] <- NA
    }
    xxnames <- xnames[fit$pivot]
    residuals <- rep.int(NA, nobs)
    residuals[good] <- z - (eta - offset)[good]
    fit$qr <- as.matrix(fit$qr)
    nr <- min(sum(good), nvars)
    if (nr < nvars) {
      Rmat <- diag(nvars)
      Rmat[1:nr, 1:nvars] <- fit$qr[1:nr, 1:nvars]
    }
    else Rmat <- fit$qr[1:nvars, 1:nvars]
    Rmat <- as.matrix(Rmat)
    Rmat[row(Rmat) > col(Rmat)] <- 0
    names(coef) <- xnames
    colnames(fit$qr) <- xxnames
    dimnames(Rmat) <- list(xxnames, xxnames)
  }
  names(residuals) <- ynames
  names(mu) <- ynames
  names(eta) <- ynames
  wt <- rep.int(0, nobs)
  wt[good] <- w^2
  names(wt) <- ynames
  names(weights) <- ynames
  names(y) <- ynames
  wtdmu <- if (intercept)
      sum(weights * y)/sum(weights)
    else linkinv(offset)
      nulldev <- sum(dev.resids(y, wtdmu, weights))
  n.ok <- nobs - sum(weights == 0)
  nulldf <- n.ok - as.integer(intercept)
  rank <- if (EMPTY) 0
  else fit$rank
  resdf <- n.ok - rank
  resdf <- n.ok
  aic.model <- aic(y, n, mu, weights, dev) + 2 * rank
  list(coefficients = coef, residuals = residuals, fitted.values = mu,
    effects = if (!EMPTY) fit$effects, R = if (!EMPTY) Rmat,
    rank = rank, qr = if (!EMPTY) structure(fit[c("qr", "rank", "qraux", "pivot", "tol")], 
    class = "qr"), family = family,
    linear.predictors = eta, deviance = dev, aic = aic.model,
    null.deviance = nulldev, iter = iter, weights = wt, prior.weights = weights,
    df.residual = resdf, df.null = nulldf, y = y, converged = conv,
    boundary = boundary, prior.mean = prior.mean, prior.scale = prior.scale,
    prior.df = prior.df, prior.sd = prior.sd, dispersion = dispersion)
}


setMethod("print", signature(x = "bayesglm"), 
    function(x, digits=2) display(object=x, digits=digits))

setMethod("show", signature(object = "bayesglm"), 
    function(object) display(object, digits=2))
