#' Single Nearest Neighborhood Matching
#'
#' Function for processing matching with propensity score
#'
#' Function for matching each treatment unit in turn the control unit
#' (not previously chosen) with the closest propensity score
#' 
#' @param z vector of indicators for treatment or control.
#' @param score vector of the propensity scores in the same order as z.
#' @param replace whether the control units could be reused for
#' matching, default is \code{FALSE}.
#' @return The function returns a vector of indices that the
#' corresponding unit is matched to. 0 means matched to nothing.
#' @references Andrew Gelman and Jennifer Hill. (2006). \emph{Data
#' Analysis Using Regression and Multilevel/Hierarchical
#' Models}. Cambridge University Press.
#' @author
#'   Jeniffer Hill \email{jh1030@@columbia.edu};
#'   Yu-Sung Su \email{suyusung@@tsinghua.edu.cn}
#' @seealso \code{\link{balance}}
#' @keywords models methods
#' @export
#' @examples
#' # matching first
#' data(lalonde)
#' attach(lalonde)
#' fit <- glm(treat ~ re74 + re75 + age + factor(educ) + 
#'             black + hisp + married + nodegr + u74 + u75, 
#'             family=binomial(link="logit"))
#' pscores <- predict(fit, type="response")
#' matches <- matching(z=lalonde$treat, score=pscores)
#' matched <- lalonde[matches$matched,]
#' 
#' # balance check!
#' b.stats <- balance(lalonde, matched, fit)
#' print(b.stats)
#' plot(b.stats)
matching <- function(z, score, replace=FALSE){
  # argument z is the vector of indicators for treatment or control #
  # argument score is the vector of the propensity scores in the    #
  # same order as z                                                 #
  # the function returns a vector of indices that the corresponding #
  # unit is matched to. 0 means matched to nothing.                 #
  #                                                                 #
  # now also returns a number for each pair making it easier to     #
  # later recover those pairs

  if (replace){
    nt <- sum(z)
    nc <- length(z) - nt
    cnts <- rep(0, nc)
    scorec <- score[z == 0]
    scoret <- score[z == 1]
    indc <- NULL
    nearest <- rep(NA, nt)
    ind.mt <- matrix(0, nc, nt)
    ind.t <- (1:(nt + nc))[z == 1]
    for(j in 1:nt) {
      near <- (1:nc)[abs(scoret[j] - scorec) == min(abs(scoret[j] - scorec))]
      if(length(near) == 1) {
        nearest[j] <- near
        indc <- c(indc, near)
      }
      else {
        nearest[j] <- near[sample(1:length(near), 1, replace = F)]
        indc <- c(indc, nearest[j])
      }
      cnts[nearest[j]] <- cnts[nearest[j]] + 1
      ind.mt[nearest[j], cnts[nearest[j]]] <- ind.t[j]
    }
    #
    ind.mt <- ind.mt[ind.mt[, 1] != 0, 1:max(cnts)]
    # now create list of indicators to pull off appropriate dataset
    ind <- numeric(nt + sum(cnts))
    # first get treat indicators
    ind[1:nt] <- (1:(nt + nc))[z == 1]
    #now the control indicators
    tmp <- (1:(nt + nc))[z == 0]
    ind[(nt + 1):length(ind)] <- tmp[indc]
    #
    out <- list(matched = unique(ind), pairs = matrix(ind, length(ind)/2, 2), ind.mt = ind.mt, cnts = cnts)
  }
  if (!replace){
    n <- length(score)
    matched <- rep(0, n)
    pairs <- rep(0, n)
    b <- (sum(z) < n/2) * 1
    tally <- 0
    for (i in (1:n)[z == b]) {
        available <- (1:n)[(z != b) & (matched == 0)]
        j <- available[order(abs(score[available] - score[i]))[1]]
        matched[i] <- j
        matched[j] <- i
        tally <- tally + 1
        pairs[c(i, j)] <- tally
    }
    out <- cbind.data.frame(matched = matched, pairs = pairs)
  }
  return(out)
}


#pscores.fun <- function(treat=Z, outs=Y, covs=X){
#  #
#  N <- nrow(covs)
#  nouts <- 1 
#  ncovs <- ncol(covs)
#  #
#  # first set up places to store results
#  res <- matrix(0,nouts,2)
#  bal <- matrix(0,ncovs,2)
#  #
#  # estimate p-scores
#  dat <- cbind.data.frame(treat=treat,covs)
#  mod <- glm(dat,family=binomial(link="logit"))
#  qx <- predict(mod, type="response")#mod$linear 
#  #
#  ### Now Matching With Replacement
#  matchout <- matching(z=treat, score=qx, replace=TRUE)
#  #
#  ### and treatment effect estimation with robust s.e.'s
#  wts <- rep(1, N)
#  wts[treat == 0] <- matchout$cnts
#  res <- .wls.all2(cbind(rep(1, sum(wts > 0)), treat[wts > 0],covs[wts > 0,  ]), wts[wts > 0], outs[wts > 0], treat[wts > 0])
#  c(res[3],sqrt(res[2]))
#}
