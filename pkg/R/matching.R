matching <- function(z, score){
  # argument z is the vector of indicators for treatment or control #
  # argument score is the vector of the propensity scores in the    #
  # same order as z                                                 #
  # the function returns a vector of indices that the corresponding #
  # unit is matched to. 0 means matched to nothing.                 #
  #                                                                 #
  # now also returns a number for each pair making it easier to     #
  # later recover those pairs
#  if 
  n <- length(score)
  matched <- rep(NA, n)
  pairs <- rep(NA, n)
  b <- (sum(z) < n/2.) * 1
  tally <- 0
  for(i in (1:n)[z == b]) {
      available <- (1:n)[(z != b) & (matched == 0)]
      j <- available[order(abs(score[available] - score[i]))[1.]]
      matched[i] <- j
      matched[j] <- i
      tally <- tally + 1.
      pairs[c(i, j)] <- tally
  }
  out <- cbind.data.frame(matched = matched, pairs = pairs)
  out
}
