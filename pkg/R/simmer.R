# the permutation matrix used during optimization, produced by cholmod.
# fill-reduces
getCholmodPerm <- function(model) {
  as(model@L@perm + 1, "pMatrix");
}

# permutes the block diagonal form into what lmer uses. See
# the discussion around getBlockCov
getRanefPerm <- function(model)
{
  numFactors <- model@dims[["nt"]];
  factorDims <- sapply(model@ST, nrow);
  numRepsPerBlock <- (model@Gp[-1] - model@Gp[-length(model@Gp)]) / factorDims;
  
  indices <- rep(0, sum(factorDims * numRepsPerBlock));

  index <- 1;
  offset <- 0;
  for (k in 1:numFactors) {
    for (l in 1:factorDims[k]) {
      for (j in 1:numRepsPerBlock[k]) {
        indices[index] <- offset + l + (j - 1) * factorDims[k];
        index <- index + 1;
      }
    }
    offset <- offset + numRepsPerBlock[k] * factorDims[k];
  }
  return(as(indices, "pMatrix"));
}

stToCov <- function(ST) {
  dimension <- nrow(ST);
  T <- ST;
  diag(T) <- rep(1, dimension);
  S <- diag(diag(ST), dimension);
  return(tcrossprod(T %*% S));
}

stToChol <- function(ST) {
  dimension <- nrow(ST);
  T <- ST;
  diag(T) <- rep(1, dimension);
  S <- diag(diag(ST), dimension);
  return(T %*% S);
}

# returns a covar (or decomp) in the form one might expect:
# S1.11 S1.12
# S1.21 S1.22
#             S1.11 S1.12
#             S1.21 S1.22
#                         S2.11 S2.12
#                         S2.21 S2.22
# lmer internally shuffles this around.
getBlockCov <- function(blockList, numRepsPerBlock) {
  numFactors <- length(blockList);
  factorDims <- sapply(blockList, nrow);

  # suppose blocks can be decomps, but ignore that for now
  numNonZeroes <- sum(factorDims^2 * numRepsPerBlock);
  rowIndices <- rep(0, numNonZeroes);
  colIndices <- rep(0, numNonZeroes);
  values <- rep(0, numNonZeroes);

  sparseIndex <- 1;
  upperLeftIndex <- 1;
  for (i in 1:numFactors) {
    numParams <- factorDims[i];
    numValues <- numParams^2;

    for (j in 1:numRepsPerBlock[i]) {
      sparseRange <- sparseIndex + 1:numValues - 1;
      upperLeftRange <- upperLeftIndex + 1:numParams - 1;

      rowIndices[sparseRange] <- rep(upperLeftRange, numParams);
      colIndices[sparseRange] <- rep(upperLeftRange, rep(numParams, numParams));
      values[sparseRange] <- as.vector(blockList[[i]]);

      sparseIndex <- sparseIndex + numValues;
      upperLeftIndex <- upperLeftIndex + numParams;
    }
  }

  return(sparseMatrix(rowIndices, colIndices, x = values));
}

# block replicates a "Sigma" for the modeled coefficients
getRanefVcov <- function(model) {
  Sigmas <- lapply(model@ST, stToCov);
  factorDims <- sapply(model@ST, nrow);
  
  numRepsPerBlock <- (model@Gp[-1] - model@Gp[-length(model@Gp)]) / factorDims;

  return(getBlockCov(Sigmas, numRepsPerBlock));
}

# block replicates a lower triangular cholesky for the modeled coefficients
getRanefChol <- function(model) {
  Lambdas <- lapply(model@ST, stToChol);
  factorDims <- sapply(model@ST, nrow);
  
  numRepsPerBlock <- (model@Gp[-1] - model@Gp[-length(model@Gp)]) / factorDims;

  return(getBlockCov(Lambdas, numRepsPerBlock));
}

# information is conditional on hyperparameters
# information is of [ranef, fixef]
getInverseInformation <- function(model) {
  # Xt.sp = sparse part of design matrix, transposed
  # X.dn  = dense part
  # transformations will be *not* applied, aside from weighting
  # observations
  # Xt.sp = L' Z' dmu/deta W^0.5
  # X.dn  = W^0.5 dmu/deta X
  #
  # dmu/deta is 1 if linear, is what it is for glm (depends on link)
  # W is 1 if linear, 1 / var(y) if glm

  Xt.sp <- model@A;
  X.dn  <- model@X;
  if (length(model@sqrtXWt) > 0) {
    if (all(dim(model@Cm)) > 0) {
      Xt.sp <- model@Cm;
    } else {
      Xt.sp@x <- model@Cx;
    }
    X.dn <- diag(as.vector(model@sqrtXWt)) %*% X.dn;
  }

  # at this point, we have the components of the augmented design matrix that
  # we need to take the crossproduct of, and then invert.

  # block-wise matrix inversion;
  # http://en.wikipedia.org/wiki/Invertible_matrix#Blockwise_inversion
  # A, B, C, and D correspond to there, although C = B' so we don't need
  # to create it
  #
  # To figure out where this comes from, you can check the docs for
  # blmer, there should be some calculations there.
  A <- tcrossprod(Xt.sp) + Diagonal(nrow(Xt.sp));
  B <- Xt.sp %*% X.dn;
  D <- crossprod(X.dn);

  A.inv <- solve(A);

  # X.i here are blocks of the inverse
  temp <- crossprod(B, A.inv); # CA^-1, comes up a lot
  D.i <- solve(D - temp %*% B);
  C.i <- -1 * D.i %*% temp;
  A.i <- A.inv - crossprod(temp, C.i);

  # since everything above is "spherical", we have to bring it to the
  # appropriate scale/rotation
  P.ranef <- getRanefPerm(model);
  Lambda <- P.ranef %*% getRanefChol(model) %*% t(P.ranef);
  
  A.i <- Lambda %*% tcrossprod(A.i, Lambda);
  C.i <- tcrossprod(C.i, Lambda);

  numRanef <- model@dims[["q"]];
  numFixef <- model@dims[["p"]];
  result <- matrix(NA, numRanef + numFixef, numRanef + numFixef);

  ranefRange <- 1:numRanef;
  fixefRange <- numRanef + 1:numFixef;
  
  result[ranefRange, ranefRange] <- as.matrix(A.i);
  result[fixefRange, ranefRange] <- as.matrix(C.i);
  result[ranefRange, fixefRange] <- t(result[fixefRange, ranefRange]);
  result[fixefRange, fixefRange] <- as.matrix(D.i);
  
  return(result);
}

# assumes p(sigma^2) propto sigma^-2
sampleSigma <- function(model) {
  return(sqrt(1 / rgamma(1, (model@dims[["n"]] - model@dims[["p"]]) / 2,
                         model@deviance[["pwrss"]] / 2)));
}

setMethod("sim", signature(object = "mer"),
          function(object, n.sims=100)
{
  # create some object matrix and util functions#########
  numRanef <- object@dims[["q"]]
  numFixef <- object@dims[["p"]]
  
  # beta.unmodeled
  fixefLabel <- names(object@fixef)
  beta.unmodeled <- NULL#array(NA, c(n.sims, numFixef))
  beta.unmodeled[[1]] <- matrix(NA, n.sims, numFixef)
  names(beta.unmodeled) <- "fixef"
  colnames(beta.unmodeled[[1]]) <- fixefLabel
  
  # beta.bygroup
  nt <- attr(object, "dim")["nt"]
  cn <- lapply(object@ST, colnames)
  ranefStructure <- lme4:::plist(lme4:::reinds(object@Gp), cn)
  wt <- lme4:::whichterms(object)
  ranefLabel <- lme4:::plist(wt, object@flist)
  grpRanefLabelIdx <- sapply(wt, length)
  grpRanefLabel <- rep(names(ranefLabel), grpRanefLabelIdx)
  
  beta.bygroup <- NULL
  for(m in 1:nt){
    rowRanefLabel <- levels(ranefLabel[[grpRanefLabel[m]]][[2]])
    colRanefLabel <- ranefStructure[[m]][[2]]
    K <- length(ranefStructure[[m]][[2]])
    J <- length(ranefStructure[[m]][[1]])/K
    beta.bygroup[[m]] <- array(NA, c(n.sims, J, K))
    dimnames(beta.bygroup[[m]]) <- c(list(NULL), list(rowRanefLabel, colRanefLabel))
  }
  names(beta.bygroup) <- grpRanefLabel
  ###################################
  
  isLinearMixedModel <- all(c(dim(object@V) == 0, length(object@muEta) == 0));

  # sigma
  simulatedSD <- ifelse(isLinearMixedModel, rep(NA, n.sims), NA);
  
  effectsMean <- c(object@ranef, object@fixef);
  # note that if anyone was counting, we could do what lmer does internally
  # and solve only the lower-triangular form against the entire
  # matrix of simulations at once. Might be faster, is essentially what
  # mvrnorm has to do internally.
  effectsCovariance <- getInverseInformation(object);

  fixefRange <- numRanef + 1:numFixef;
  
  for (i in 1:n.sims) {
    if (isLinearMixedModel) {
      simulatedSD[i] <- sampleSigma(object);
        # ifelse(object@dims[["REML"]] == 1, object@deviance[["sigmaREML"]], object@deviance[["sigmaML"]]);
      simulatedEffects <- mvrnorm(1, effectsMean, (simulatedSD[i] ^ 2) * effectsCovariance)
    } else {
      simulatedEffects <- mvrnorm(1, effectsMean, effectsCovariance);
    }
      
    beta.unmodeled[[1]][i,] <- simulatedEffects[fixefRange]
    for (m in 1:nt){
      dim.ranef <- dim(beta.bygroup[[m]])
      idx <- ranefStructure[[m]][[1]]
      beta.bygroup[[m]][i,,] <- matrix(simulatedEffects[idx], dim.ranef[2], dim.ranef[3])
    }
  }
  ###############
  # reorganize ranef
  reGrpRanef <- function(el, ...){
    ans <- do.call(abind, beta.bygroup[el[[1]]])
  } 
  beta.bygroup <- lapply(ranefLabel, reGrpRanef, beta.bygroup)
  rr <- ranef(object)
  n.grp <- length(rr)
  for(m in 1:n.grp){
    dimnames(beta.bygroup[[m]]) <- c(list(NULL), dimnames(rr[[m]]))
  }
  #####################

  ans <- new("sim.mer", 
             "fixef" = beta.unmodeled$fixef,
             "ranef" = beta.bygroup,
             "sigma" = simulatedSD)
  return(ans)
})
