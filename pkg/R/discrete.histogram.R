#' Histogram for Discrete Distributions
#'
#' Creates a prettier histogram for discrete distributions
#'
#' This function displays a histogram for discrete probability
#' distributions.
#'
#' @param x The vector of x's.
#' @param prob The probabilities for the x's.
#' @param prob2 A second vector of probabilities of the x's.
#' @param prob3 A third vector of probabilities of the x's.
#' @param xlab Label for the x axis.
#' @param xaxs.label Label for the x's.
#' @param yaxs.label Label for the y axis.
#' @param bar.width Width of the bars.
#' @param freq If TRUE, shows a frequency histogram as opposed to
#'   probability.
#' @param prob.col The color of the first set of histogram bars.
#' @param prob2.col The color of the second set of histogram bars.
#' @param prob3.col The color of the third set of histogram bars.
#' @param ... Additional arguments passed to function \code{plot}.
#' 
#' @keywords dplot
#' @examples
#' a <- c(3,4,0,0,5,1,1,1,1,0)
#' discrete.histogram (a)
#' 
#' x <- c(0,1,3,4,5)
#' p <- c(.3,.4,.1,.1,.1)
#' discrete.histogram (x,p)
#' 
#' x <- c(0,1,3,4,5)
#' y <- c(3,4,1,1,1)
#' discrete.histogram (x,y)
discrete.histogram <- function (x, prob, prob2 = NULL, prob3 = NULL, xlab = "x", xaxs.label = NULL,
    yaxs.label = NULL, bar.width = NULL, freq = FALSE, prob.col = "blue",
    prob2.col = "red", prob3.col = "gray", ...)
{
    if (!missing(x) && missing(prob)) {
        prob <- table(x)
        x <- sort(unique(x))
    }
    if (length(x) != length(prob)) {
        stop("Length of 'x' must be the same as the length of 'prob'")
    }
    if (!freq) {
        prob <- prob/sum(prob)
        prob2 <- prob2/sum(prob2)
        prob3 <- prob3/sum(prob3)
        ylab <- "Probability"
    }
    else {
        ylab <- "Count"
    }
    if (is.numeric(x)) {
        x.values <- sort(unique(x))
        n.x.values <- length(x.values)
        if (is.null(bar.width)) {
            gaps <- x.values[2:n.x.values] - x.values[1:(n.x.values -
                1)]
            bar.width <- min(gaps) * 0.2
        }
        par(mar = c(3, 3, 4, 1), mgp = c(1.7, 0.5, 0), tck = -0.01)
        plot(range(x) + c(-2, 2) * bar.width, c(0, max(prob,
            prob2, prob3)), xlab = xlab, ylab = ylab, xaxs = "i",
            xaxt = "n", yaxs = "i", yaxt = ifelse(is.null(yaxs.label),
                "s", "n"), bty = "l", type = "n", ...)
        if (is.null(xaxs.label)) {
            axis(1, x.values)
        }
        else {
            axis(1, xaxs.label[[1]], xaxs.label[[2]])
        }
    }
    else {
        x.values <- unique(x)
        n.x.values <- length(x.values)
        if (is.null(bar.width)) {
            bar.width <- 0.2
        }
        par(mar = c(3, 3, 4, 1), mgp = c(1.7, 0.5, 0), tck = -0.01)
        plot(c(1, n.x.values) + c(-2, 2) * bar.width, c(0, max(prob,
            prob2, prob3)), xlab = xlab, ylab = ylab, xaxs = "i",
            xaxt = "n", yaxs = "i", yaxt = ifelse(is.null(yaxs.label),
                "s", "n"), bty = "l", type = "n", ...)
        if (is.null(xaxs.label)) {
            axis(1, 1:n.x.values, x.values)
        }
        else {
            axis(1, xaxs.label[[1]], xaxs.label[[2]])
        }
        x <- 1:length(x)
    }
    if (!is.null(yaxs.label)) {
        axis(2, yaxs.label[[1]], yaxs.label[[2]])
    }
    offset <- rep(0, 3)
    if (length(prob2) != 0 & length(prob3) != 0) {
        offset[1] <- -bar.width
        offset[2] <- 0
        offset[3] <- bar.width
    }
    if (length(prob2) > 0 & length(prob3) == 0) {
        offset[1] <- -bar.width/2
        offset[2] <- bar.width/2
        offset[3] <- 0
    }
    for (i in 1:length(x)) {
        polygon(x[i] + c(-1, -1, 1, 1) * bar.width/2 + offset[1],
            c(0, prob[i], prob[i], 0), border = prob.col, col = prob.col)
        if (!is.null(prob2)) {
            polygon(x[i] + c(-1, -1, 1, 1) * bar.width/2 + offset[2],
                c(0, prob2[i], prob2[i], 0), border = prob2.col,
                col = prob2.col)
        }
        if (!is.null(prob3)) {
            polygon(x[i] + c(-1, -1, 1, 1) * bar.width/2 + offset[3],
                c(0, prob3[i], prob3[i], 0), border = prob3.col,
                col = prob3.col)
        }
    }
}

#' @rdname discrete.histogram
#' @export
discrete.hist <- discrete.histogram
