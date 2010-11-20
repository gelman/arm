getAugmentedAndRotatedDesignMatrix <- function(model) {
  numModeledParameters <- model@dims[['q']];
  numUnmodeledParameters <- model@dims[['p']];
  numObservations <- model@dims[['n']];
  
  numParameters <- numUnmodeledParameters + numModeledParameters;
  numAugmentedRows <- numModeledParameters + numObservations;
  
  # set up design matrix:
  # [ ZLambda X ]
  # [ I       0 ]
  designMatrix <- Matrix(0, numAugmentedRows, numParameters, sparse=TRUE);
  if (numModeledParameters > 0) {
    designMatrix[1:numObservations, 1:numModeledParameters] <- t(model@A);
    designMatrix[(1 + numObservations):numAugmentedRows, 1:numModeledParameters] <- diag(1, numModeledParameters);
  }
  if (numUnmodeledParameters > 0) {
    designMatrix[1:numObservations, (numModeledParameters + 1):numParameters] <- model@X;
  }

  return(designMatrix);
}

getAugmentedDesignMatrix <- function(model) {
  numModeledParameters <- model@dims[['q']];
  numUnmodeledParameters <- model@dims[['p']];
  numObservations <- model@dims[['n']];
  
  numParameters <- numUnmodeledParameters + numModeledParameters;
  numAugmentedRows <- numModeledParameters + numObservations;
  
  # set up design matrix:
  # [ Z          X ]
  # [ Lambda^-1  0 ]
  designMatrix <- Matrix(0, numAugmentedRows, numParameters, sparse=TRUE);
  if (numModeledParameters > 0) {
    designMatrix[1:numObservations, 1:numModeledParameters] <- t(model@Zt);
    designMatrix[(1 + numObservations):numAugmentedRows, 1:numModeledParameters] <-
      solve(getSphericalCovarianceLeftFactor(model));
  }
  if (numUnmodeledParameters > 0) {
    designMatrix[1:numObservations, (numModeledParameters + 1):numParameters] <- model@X;
  }

  return(designMatrix);
}

getAugmentedResponse <- function(model) {
  numModeledParameters <- model@dims[['q']];
  numObservations <- model@dims[['n']];

  numAugmentedRows <- numModeledParameters + numObservations;

  response <- Matrix(0, numAugmentedRows, sparse=TRUE);
  response[1:numObservations] <- model@y;

  return(response);
}

getNumModeledParametersPerLevel <- function(model) {
  numLevels <- model@dims[["nt"]];

  return(sapply(model@ST, nrow))
}

getNumGroupsPerLevel <- function(model) {
  numLevels <- model@dims[["nt"]];
  return (sapply(1:numLevels, function(i) {
    return ((model@Gp[i + 1] - model@Gp[i]) / nrow(model@ST[[i]]));
  }));
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

  numLevels <- model@dims[["nt"]];
  numModeledParametersPerLevel <- getNumModeledParametersPerLevel(model);
  numGroupsPerLevel <- getNumGroupsPerLevel(model);
  
  numModeledParameters <- model@dims[["q"]];

  covarianceMatrix <- Matrix(0, numModeledParameters, numModeledParameters, sparse=TRUE);
  
  covarianceIndex <- 1;
  for (k in 1:numLevels) {
    # store L = T %*% S
    triangularMatrix <- model@ST[[k]];
    diag(triangularMatrix) <- rep(1, numModeledParametersPerLevel[k]);
    
    diagonalMatrix <- diag(diag(model@ST[[k]]), numModeledParametersPerLevel[k]);

    leftFactor <- triangularMatrix %*% diagonalMatrix;

    # copy in L into covariance as many times as needed, one block for each group
    # at this level
    for (j in 1:numGroupsPerLevel[k]) {
      firstRow <- covarianceIndex;
      lastRow <- covarianceIndex + numModeledParametersPerLevel[k] - 1;
      
      covarianceMatrix[firstRow:lastRow, firstRow:lastRow] <- triangularMatrix %*% diagonalMatrix;
      
      covarianceIndex <- covarianceIndex + numModeledParametersPerLevel[k];
    }
  }
  
  return(covarianceMatrix);
}

# this excludes the common scale factor, which it is assumed will be
# simulated on its own
getLinearCovarianceEstimate <- function(model) {
  numObservations <- model@dims[["n"]];
  numModeledParameters <- model@dims[["q"]];
  
  numAugmentedRows <- numObservations + numModeledParameters;
  
  designMatrix <- getAugmentedDesignMatrix(model);

  return(solve(crossprod(designMatrix)));
  
  #responseCovariance <- Matrix(0, numAugmentedRows, numAugmentedRows, sparse=TRUE);
  #responseCovariance[1:numObservations, 1:numObservations] <-
  #  Diagonal(numObservations, 1) + crossprod(model@A);
  #responseCovariance <- Matrix(responseCovariance, sparse=TRUE);
  #designCrossprodInverse <- solve(crossprod(designMatrix));

  # (X'X)^-1 X' covar([y, 0]') X (X'X)^-1',
  # covar([ y ]) = [ I.n + ZL L'Z'  0 ]
  #      ([ 0 ]) = [       0        0 ]
  #return (designCrossprodInverse %*% t(designMatrix) %*% responseCovariance %*%
  #        designMatrix %*% designCrossprodInverse);
}

getPermutationMatrix <- function(model) {
  return(as(model@L@perm + 1, "pMatrix"));
}

# information is conditional on hyperparameters
# information is of [ranef, fixef]
getObservedInformation <- function(model) {
  numModeledParameters <- model@dims[["q"]];
  numUnmodeledParameters <- model@dims[["p"]];

  numParameters <- numModeledParameters + numUnmodeledParameters;
  
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
  upperLeftBlock <- tcrossprod(ranefWeightedJacobianTrans) + diag(1, numModeledParameters);
  lowerRightBlock <- crossprod(fixefWeightedJacobian);
  upperRightBlock <- ranefWeightedJacobianTrans %*% fixefWeightedJacobian;
  lowerLeftBlock <- t(upperRightBlock);

  observedInformation <- Matrix(0, numParameters, numParameters, sparse=TRUE);
  observedInformation[1:numModeledParameters, 1:numModeledParameters] <-
    upperLeftBlock;
  observedInformation[(numModeledParameters + 1):numParameters, (numModeledParameters + 1):numParameters] <-
    lowerRightBlock;
  observedInformation[1:numModeledParameters, (numModeledParameters + 1):numParameters] <-
    upperRightBlock;
  observedInformation[(numModeledParameters + 1):numParameters, 1:numModeledParameters] <-
    lowerLeftBlock;

  return(observedInformation);
}

sampleSigma <- function(model) {
  # proportional to a chi-square with d.o.f. = num observations + num ranef
  return(1 / rgamma(1, model@dims[["n"]] / 2, model@deviance[["pwrss"]] / 2));
  
  return(sqrt(model@deviance[["pwrss"]] / rchisq(1, model@dims[["n"]])));# + model@dims[["q"]])));
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
  #rotation <- Diagonal(numRanef + numFixef, 1);
  #rotation[1:numRanef, 1:numRanef] <- getSphericalCovarianceLeftFactor(object) %*% t(getPermutationMatrix(object));
  #effectsCovariance <- rotation %*% tcrossprod(effectsCovariance, rotation);

  if (isLinearMixedModel) {
    for (i in 1:n.sims) {
      simulatedSD[[1]][i] <- #sampleSigma(object);
        ifelse(object@dims[["REML"]] == 1, object@deviance[["sigmaREML"]], object@deviance[["sigmaML"]]);
      simulatedEffects <- mvrnorm(1, effectsMean, (simulatedSD[[1]][i] ^ 2) * effectsCovariance);
      beta.unmodeled[[1]][i,] <- simulatedEffects[(numRanef + 1):(numRanef + numFixef)]
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
                "fixef" = beta.unmodeled$fixef,
                "ranef" = beta.bygroup,
                "sigma" = NA)

    return(ans)
  }
})
