#' Logistic and Inverse logistic functions
#'
#' Inverse-logit function, transforms continuous values to the range
#' (0, 1)
#'
#' The Inverse-logit function defined as: \eqn{logit^-1(x) =
#' e^x/(1+e^x)} transforms continuous values to the range (0, 1),
#' which is necessary, since probabilities must be between 0 and 1 and
#' maps from the linear predictor to the probabilities
#'
#' @param x A vector of continuous values.
#' @return A vector of estimated probabilities.
#' @references Andrew Gelman and Jennifer Hill. (2006). \emph{Data
#'   Analysis Using Regression and Multilevel/Hierarchical
#'   Models}. Cambridge University Press.
#' @author
#'   Andrew Gelman \email{gelman@@stat.columbia.edu},
#'   M.Grazia Pittau \email{grazia@@stat.columbia.edu}
#' @keywords models
#' @export
#' @examples
#' data(frisk)
#' n <- 100
#' x1 <- rnorm (n)
#' x2 <- rbinom (n, 1, .5)
#' b0 <- 1
#' b1 <- 1.5
#' b2 <- 2
#' Inv.logit <- invlogit(b0+b1*x1+b2*x2)
#' plot(b0+b1*x1+b2*x2, Inv.logit)
invlogit <- function (x) {
    1/(1+exp(-x))
}

#' @rdname invlogit
#' @export
logit <- function (x) {
  log(x/(1-x))
}
