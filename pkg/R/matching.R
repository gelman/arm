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
    ind.mt <- ind.mt[ind.mt[, 1], 1:max(cnts)]
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
