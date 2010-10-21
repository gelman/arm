getAugmentedDesignMatrix <- function(model) {
  numRanef <- model@dims[['q']];
  numFixef <- model@dims[['p']];
  numObservations <- model@dims[['n']];
  # set up design matrix:
  # [ ZLambda X ]
  # [ I       0 ]
  designMatrix <- Matrix(0, numObservations + numRanef, numRanef + numFixef, sparse=TRUE);
  if (numRanef > 0) {
    designMatrix[1:numObservations, 1:numRanef] <- t(model@A);
    designMatrix[(1 + numObservations):(numObservations + numRanef), 1:numRanef] <- diag(1, numRanef);
  }
  if (numFixef > 0) {
    designMatrix[1:numObservations, (numRanef + 1):(numRanef + numFixef)] <- model@X;
  }

  return(designMatrix);
}

getAugmentedResponse <- function(model) {
  numRanef <- model@dims[['q']];
  numObservations <- model@dims[['n']];

  response <- Matrix(0, numRanef + numObservations, sparse=TRUE);
  response[1:numObservations] <- model@y;

  return(response);
}

getNumColumnsPerTerm <- function(model) {
  numTerms <- length(model@ST);
  numColumns <- rep(0, numTerms);
  
  for (i in 1: numTerms) {
    numColumns[i] <- dim(model@ST[[i]])[1];
  }
  return(numColumns);
}

getNumLevelsPerTerm <- function(model) {
  numTerms <- length(model@ST);
  numLevels <- rep(0, numTerms);
  
  for (i in 1:numTerms) {
    numLevels[i] <- (model@Gp[i + 1] - model@Gp[i]) / dim(model@ST[[i]])[1];
  }
  return(numLevels);
}

getTermIndex <- function(index, numTerms, termPointers) {
  for (i in 1:numTerms){
    if (index < termPointers[i + 1]){
      return (i + 1);
    }
    else{
      return (-1);
    }
  }
}

getSphericalCovarianceLeftFactor <- function(model) {
  # is T %*% S, but we have to fill T and S from the block form
  # in which they are stored. Unfortunately, the blocks cannot
  # simply be replicated, as the random effects are not stored
  # in block form themselves (blocks are per level, not per term)

  # should be lower triangular

  numTerms <- model@dims[["nt"]];
  numRanef <- model@dims[["q"]];
  termPointers <- model@Gp;

  numLevelsPerTerm <- getNumLevelsPerTerm(model);
  numColumnsPerTerm <- getNumColumnsPerTerm(model);

  triangularMatrix <- diag(1, numRanef);
  diagonalMatrix <- Matrix(0, numRanef, numRanef, sparse=TRUE);
  covarianceMatrix <- Matrix(0, numRanef, numRanef, sparse=TRUE);
  
  covarianceIndex <- 1;
  for (i in 1:numTerms) {
    triangularMatrix <- model@ST[[i]];
    diag(triangularMatrix) <- rep(1, dim(model@ST[[i]])[1]);
    diagonalMatrix <- diag(diag(model@ST[[i]]), dim(model@ST[[i]])[1]);
    
    for (j in 1:numLevelsPerTerm[i]) {
      firstRow <- covarianceIndex;
      lastRow <- covarianceIndex + numColumnsPerTerm[i] - 1;
      
      covarianceMatrix[firstRow:lastRow, firstRow:lastRow] <- triangularMatrix %*% diagonalMatrix;
      
      covarianceIndex <- covarianceIndex + numColumnsPerTerm[i];
    }
  }
  
  return(covarianceMatrix);
}

# this excludes the common scale factor, which it is assumed will be
# simulated on its own
getLinearCovarianceEstimate <- function(model) {
  numObservations <- model@dims[["n"]];
  numRanef <- model@dims[["q"]];
  
  designMatrix <- getAugmentedDesignMatrix(model);
  responseCovariance <- Matrix(0, numObservations + numRanef, numObservations + numRanef, sparse=TRUE);
  responseCovariance[1:numObservations, 1:numObservations] <- Diagonal(numObservations, 1) + crossprod(model@A);
  responseCovariance <- Matrix(responseCovariance, sparse=TRUE)
  designCrossprodInverse <- solve(crossprod(designMatrix));

  # (X'X)^-1 X' covar([y 0]') X (X'X)^-1'
  return (designCrossprodInverse %*% t(designMatrix) %*% responseCovariance %*% designMatrix %*% designCrossprodInverse);
}

getPermutationMatrix <- function(model) {
  return(as(model@L@perm + 1, "pMatrix"));
}

# information is conditional on hyperparameters
# information is of [ranef, fixef]
getObservedInformation <- function(model) {
  numRanef <- model@dims[["q"]];
  numFixef <- model@dims[["p"]];
  
  # the matrix C, which corresponds to ST'Z'diag(d mu / d eta)W^.5,
  # is closely related to the update matrix for the random effects
  # in literature, this is U'
  #ranefWeightedJacobianTrans <- NULL;
  #sparseWeightMatrix <- NULL;
  if (length(model@sqrtXWt) > 0) {
    if (all(dim(model@Cm) > 0)) {
      ranefWeightedJacobianTrans <- model@Cm;
    } 
    else {
      # this should be obtained from model@Cx, but that requires rebuilding
      # the matrix from A. For now, it is re-computed
      sparseWeightMatrix <- Diagonal(length(as.vector(model@sqrtXWt)), as.vector(model@sqrtXWt));
      
      ranefWeightedJacobianTrans <- model@A %*% sparseWeightMatrix;
    }
  } 
  else {
    ranefWeightedJacobianTrans <- model@A;
  }

  # matrix for beta is more difficult
  # X (d mu / d eta) W ^.5
  # in literature, this is V
  #fixefWeightedJacobian <- NULL;
  if (all(dim(model@sqrtXWt) > 0)) {
    if (is.null(sparseWeightMatrix)) {
      sparseWeightMatrix <- Diagonal(length(as.vector(model@sqrtXWt)), as.vector(model@sqrtXWt));
    }
    fixefWeightedJacobian <- sparseWeightMatrix %*% model@X;
  } 
  else {
    fixefWeightedJacobian <- model@X;
  }

  # if the above are U and V, the observed information is:
  # [ U'U + I   U'V ]
  # [ V'U       V'V ]
  upperLeftBlock <- tcrossprod(ranefWeightedJacobianTrans) + diag(1, numRanef);
  lowerRightBlock <- crossprod(fixefWeightedJacobian);
  upperRightBlock <- ranefWeightedJacobianTrans %*% fixefWeightedJacobian;
  lowerLeftBlock <- t(upperRightBlock);

  observedInformation <- Matrix(0, numRanef + numFixef, numRanef + numFixef, sparse=TRUE);
  observedInformation[1:numRanef, 1:numRanef] <- upperLeftBlock;
  observedInformation[(numRanef + 1):(numRanef + numFixef), (numRanef + 1):(numRanef + numFixef)] <- lowerRightBlock;
  observedInformation[1:numRanef, (numRanef + 1):(numRanef + numFixef)] <- upperRightBlock;
  observedInformation[(numRanef + 1):(numRanef + numFixef), 1:numRanef] <- lowerLeftBlock;

  return(observedInformation);
}

sampleSigma <- function(model) {
  # proportional to a chi-square with d.o.f. = num observations + num ranef
  return(sqrt(model@deviance[["pwrss"]] / rchisq(1, model@dims[["n"]] + model@dims[["q"]])));
}

setMethod("sim", signature(object = "mer"),
          function(object, n.sims=100)
{
  # create some object matrix and util functions#########
  numRanef <- object@dims[["q"]]
  numFixef <- object@dims[["p"]]
  # sigma
  simulatedSD <- NULL
  simulatedSD[[1]] <- rep(NA, n.sims)
  names(simulatedSD) <- "sigma"
  
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
  
  effectsMean <- c(object@ranef, object@fixef);
  effectsCovariance <- NULL;
  if (isLinearMixedModel) {
    effectsCovariance <- getLinearCovarianceEstimate(object);
  } else {
    effectsCovariance <- solve(getObservedInformation(object));
  }

  # the random effects in the model (b) differ from the values simulated (u in the model):
  #   b = Lambda %*% t(P) %*% u
  #
  # where Lambda is the left factor of Sigma and P is the permutation matrix calculated
  # when taking the left factor of (Lambda Z Z' Lambda' + I)
  rotation <- Diagonal(numRanef + numFixef, 1);
  rotation[1:numRanef, 1:numRanef] <- getSphericalCovarianceLeftFactor(object) %*% t(getPermutationMatrix(object));
  
  effectsCovariance <- rotation %*% effectsCovariance %*% t(rotation);
  
  #simulatedEffects <- array(0, c(n.sims, numRanef + numFixef));

  # I am not yet clear on R's ability to optimize loops, so I've moved the conditional outside
  if (isLinearMixedModel) {
    for (i in 1:n.sims) {
      simulatedSD[[1]][i] <- sampleSigma(object);
      simulatedEffects <- mvrnorm(1, effectsMean, (simulatedSD[[1]][i] ^ 2) * effectsCovariance);
      beta.unmodeled[[1]][i,] <- simulatedEffects[(numRanef + 1):(numRanef + numFixef)]
      for(m in 1:nt){
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
    browser()
    ans <- new("sim.mer", 
               "fixef" = beta.unmodeled$fixef,
               "ranef" = beta.bygroup,
               "sigma" = simulatedSD$sigma)
    return(ans)
  } 
  else {
    for (i in 1:n.sims) {
      simulatedEffects <- mvrnorm(1, effectsMean, effectsCovariance)
      beta.unmodeled[[1]][i,] <- simulatedEffects[(numRanef + 1):(numRanef + numFixef)]
      for(m in 1:nt){
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
                beta.unmodeled,
                beta.bygroup,
                "sigma" = NA)

    return(ans)
  }
})
