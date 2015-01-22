#' Formating the Rounding of Numbers
#'
#' @description
#'   \code{fround} rounds the values in its first argument to the
#'   specified number of decimal places with surrounding quotes.
#' 
#'   \code{pfround} rounds the values in its first argument to the
#'   specified number of decimal places without surrounding quotes.
#'
#' @param x a numeric vector.
#' @param digits integer indicating the precision to be used.
#' @author
#'   Andrew Gelman \email{gelman@@stat.columbia.edu};
#'   Yu-Sung Su \email{suyusung@@tsinghua.edu.cn}
#' @seealso \code{\link{round}}
#' @keywords manip print
#' @export
#' @examples
#' x <- rnorm(1)
#' fround(x, digits=2)
#' pfround(x, digits=2)
fround <- function (x, digits) {
    format (round (x, digits), nsmall=digits)
}

#' @rdname fround
#' @export  
pfround <- function (x, digits) {
    print (fround (x, digits), quote=FALSE)
}
 
