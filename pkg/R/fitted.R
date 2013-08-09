## the plan here is to shuffle the ranefs back into the way a mer object
## stores them so that a simple X * beta + Z * theta op does the trick
#setMethod("fitted", signature(object = "sim.mer"),
#          function(object, regression)
#{
#  if (missing(regression) || is.null(regression)) stop("fitted for sim.mer requires original merPred object as well.");
#  if (!inherits(regression, "mer")) stop("regression argument for fitted on sim.mer does not inherit from class 'mer'");
#  sims <- object;
#  numSimulations <- dim(sims@fixef)[1];
#  numRanef <- regression@dims[["q"]];
#  numLevels <- regression@dims[["nt"]];
#  
#  ranefNames <- lapply(regression@ST, colnames);
#  # ranefStructure <=> list w/one element per level, each is list with first element
#  # indicies into joint covar, second with name of what is varying
#  ranefStructure <- lme4:::plist(lme4:::reinds(regression@Gp), ranefNames);
#  wt <- lme4:::whichterms(regression);
#  ranefLabel <- lme4:::plist(wt, regression@flist);
#
#  ranef <- matrix(0, numRanef, numSimulations);
#  
#  levelIndex <- 1;
#  numUniqueLevels <- length(ranefLabel);
#  
#  for (uniqueLevelNum in 1:numUniqueLevels) {
#    numCombinedLevels <- length(ranefLabel[[uniqueLevelNum]][[1]]);
#
#    simIndex <- 0;
#    for (subLevelNum in 1:numCombinedLevels) {
#      ranefIndices <- ranefStructure[[levelIndex]][[1]];
#      simIndices <- 1:length(ranefStructure[[levelIndex]][[2]]) + simIndex;
#
#      numCoefs <- dim(sims@ranef[[uniqueLevelNum]])[2];
#      for (i in 1:length(simIndices)) {
#        ranefSubset <- 1:numCoefs + (i - 1) * numCoefs;
#        
#        ranef[ranefIndices[ranefSubset],] <- t(sims@ranef[[uniqueLevelNum]][,,simIndices[i]]);
#      }
#      
#      simIndex <- simIndex + length(simIndices);
#      levelIndex <- levelIndex + 1;
#    }
#  }
#
#  linkType <- regression@dims[["lTyp"]];
#  invLinks <- c(function(x) { e.x <- exp(x); e.x / (1 + e.x); }, # logit
#                function(x) { pnorm(x); },                       # probit
#                function(x) { pcauchy(x); },                    # cauchit
#                function(x) { 1 - exp(-exp(x)); },               # comp log-log
#                function(x) { x; },                              # identity
#                function(x) { exp(x); },                         # log
#                function(x) { x^2; },                            # sqrt
#                function(x) { x^-0.5; },                         # 1 / mu^2
#                function(x) { 1 / x; });                         # inverse
#  invLink <- invLinks[[linkType]];
#
#  result <- tcrossprod(regression@X, sims@fixef) + crossprod(regression@Zt, ranef);
#  return(invLink(result));
#});

# the plan here is to shuffle the ranefs back into the way a mer object
# stores them so that a simple X * beta + Z * theta op does the trick
setMethod("fitted", signature(object = "sim.merMod"),
          function(object, regression)
{
  if (missing(regression) || is.null(regression)) stop("fitted for sim.mer requires original merPred object as well.");
  if (!inherits(regression, "merMod")) stop("regression argument for fitted on sim.mer does not inherit from class 'merMod'");
  sims <- object;
  numSimulations <- dim(sims@fixef)[1];
  devcomp <- getME(regression, "devcomp");
  dims <- devcomp$dims;
  
  numRanef  <- dims[["q"]];
  numLevels <- dims[["reTrms"]];

  simulatedRanef <- matrix(0, numRanef, numSimulations);

  index <- 0;
  for (i in 1:length(sims@ranef)) {
    levelSims <- sims@ranef[[i]];
    numCoefficientsPerLevel <- dim(levelSims)[2];
    numGroupsPerLevel <- dim(levelSims)[3];
    for (j in 1:numCoefficientsPerLevel) {
      ranefRange <- index + 1:numGroupsPerLevel;
      index <- index + numGroupsPerLevel;
      
      simulatedRanef[ranefRange,] <- t(levelSims[,j,]);
    }
  }

  X <- getME(regression, "X");
  Zt <- getME(regression, "Zt");
  W <- Diagonal(nrow(X), regression@resp$sqrtXwt);
  X <- W %*% X;
  Zt <- Zt %*% W;
  
  result <- tcrossprod(X, sims@fixef) + crossprod(Zt, simulatedRanef);
  
  if (devcomp$dims[["GLMM"]] == 0L) return(result);

  return(regression@resp$family$linkinv(result));
});
