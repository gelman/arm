setMethod("fitted", signature(object = "sim.mer"),
          function(object, regression)
{
  if (missing(regression) || is.null(regression)) stop("fitted for sim.mer requires original mer object as well.");
  if (!inherits(regression, "mer")) stop("regression argument for fitted on sim.mer does not inherit from class 'mer'");
  sims <- object;
  numSimulations <- dim(sims@fixef)[1];
  numRanef <- regression@dims[["q"]];
  numLevels <- regression@dims[["nt"]];
  
  ranefNames <- lapply(regression@ST, colnames);
  # ranefStructure <=> list w/one element per level, each is list with first element
  # indicies into joint covar, second with name of what is varying
  ranefStructure <- lme4:::plist(lme4:::reinds(regression@Gp), ranefNames);
  wt <- lme4:::whichterms(regression);
  ranefLabel <- lme4:::plist(wt, regression@flist);

  ranef <- matrix(0, numRanef, numSimulations);
  
  levelIndex <- 1;
  numUniqueLevels <- length(ranefLabel);
  
  for (uniqueLevelNum in 1:numUniqueLevels) {
    numCombinedLevels <- length(ranefLabel[[uniqueLevelNum]][[1]]);

    simIndex <- 0;
    for (subLevelNum in 1:numCombinedLevels) {
      ranefIndices <- ranefStructure[[levelIndex]][[1]];
      simIndices <- 1:length(ranefStructure[[levelIndex]][[2]]) + simIndex;

      ranef[ranefIndices,] <- sims@ranef[[uniqueLevelNum]][,,simIndices];
      
      simIndex <- simIndex + length(simIndices);
      levelIndex <- levelIndex + 1;
    }
  }

  linkType <- regression@dims[["lTyp"]];
  invLinks <- c(function(x) { e.x <- exp(x); e.x / (1 + e.x); }, # logit
                function(x) { pnorm(x); },                       # probit
                function(x) { pcauchy(x); },                    # cauchit
                function(x) { 1 - exp(-exp(x)); },               # comp log-log
                function(x) { x; },                              # identity
                function(x) { exp(x); },                         # log
                function(x) { x^2; },                            # sqrt
                function(x) { x^-0.5; },                         # 1 / mu^2
                function(x) { 1 / x; });                         # inverse
  invLink <- invLinks[[linkType]];

  result <- tcrossprod(regression@X, sims@fixef) + crossprod(regression@Zt, ranef);
  return(invLink(result));
});
