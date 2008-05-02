matching <- function(z, score){

    n <- length(score)
    matched <- rep(0., n)
    pairs <- rep(0., n)
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
