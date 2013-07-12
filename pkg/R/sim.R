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
    dimnames(beta) <- list (NULL, names(beta.hat))
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
