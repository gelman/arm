#' Contrast Matrices
#'
#' Return a matrix of contrasts used in \code{\link{bayesglm}}.
#'
#' These functions are adapted from \code{contr.treatment} and
#' \code{contr.poly} in \code{\link{stats}} package.  The purpose for
#' these functions are to keep the baseline levels of categorical
#' variables and thus to suit the use of \code{\link{bayesglm}}.
#'
#' \code{contr.bayes.unordered} is equivalent to
#' \code{contr.treatment} whereas \code{contr.bayes.ordered} is
#' equivalent to \code{contr.poly}.
#'
#' @param n a vector of levels for a factor, or the number of levels.
#' @param base an integer specifying which group is considered the
#' baseline group. Ignored if \code{contrasts} is \code{FALSE}.
#' @param contrasts a logical indicating whether contrasts should be computed.
#' @param scores the set of values over which orthogonal polynomials
#' are to be computed.
#'
#' @author Yu-Sung Su \email{suyusung@@tsinghua.edu.cn}
#' @seealso \code{\link{C}}, \code{\link{contr.helmert}},
#' \code{\link{contr.poly}}, \code{\link{contr.sum}},
#' \code{\link{contr.treatment}}; \code{\link{glm}},
#' \code{\link{aov}}, \code{\link{lm}}, \code{\link{bayesglm}}.
#' @keywords design regression array manip
#' @name contrasts.bayes
#' @examples
#' cat.var <- rep(1:3, 5)
#' dim(contr.bayes.unordered(cat.var))
#' # 15*15 baseline level kept!
#' dim(contr.treatment(cat.var))
#' # 15*14
NULL

#' @rdname contrasts.bayes
#' @export
contr.bayes.ordered <- function ( n, scores = 1:n, contrasts = TRUE )
{
    make.poly <- function( n, scores ) {
        y   <- scores - mean( scores )
        X   <- outer( y, seq_len( n ) - 1, "^" )
        QR  <- qr( X )
        z   <- QR$qr
        z   <- z *( row( z ) == col( z ) )
        raw <- qr.qy( QR, z )
        Z   <- sweep( raw, 2, apply( raw, 2, function( x ) sqrt( sum( x^2 ) ) ), "/" )
        colnames( Z ) <- paste( "^", 1:n - 1, sep="" )
        Z
    }
    if ( is.numeric( n ) && length( n ) == 1 ) { levs <- 1:n }
    else {
        levs <- n
        n <- length( levs )
    }
    if ( n < 2 ) {
        stop( gettextf( "contrasts not defined for %d degrees of freedom", n - 1 ), domain = NA ) 
    }
    if ( n > 95 ) {
        stop( gettextf( "orthogonal polynomials cannot be represented accurately enough for %d degrees of freedom", n-1 ), domain = NA ) 
    }
    if ( length( scores ) != n ) {
        stop( "'scores' argument is of the wrong length" )
    }
    if ( !is.numeric( scores ) || any( duplicated( scores ) ) ) {
        stop("'scores' must all be different numbers")
    }
    contr <- make.poly( n, scores )
    if ( contrasts ) {
        dn <- colnames( contr )
        dn[2:min( 4, n )] <- c( ".L", ".Q", ".C" )[1:min( 3, n-1 )]
        colnames( contr ) <- dn
        contr[, , drop = FALSE]
    }
    else {
        contr[, 1] <- 1
        contr
    }
}

#' @rdname contrasts.bayes
#' @export
contr.bayes.unordered <- function(n, base = 1, contrasts = TRUE) {
    if( is.numeric( n ) && length( n ) == 1) {
        if( n > 1 ) { levs <- 1:n }
        else stop( "not enough degrees of freedom to define contrasts" )
    } 
    else {
        levs <- n
        n <- length( n )
    }
    contr <- array( 0, c(n, n), list( levs, levs ) )
    diag( contr ) <- 1
    if( contrasts ) {
        if( n < 2 ) { stop( gettextf( "contrasts not defined for %d degrees of freedom", n - 1 ), domain = NA ) }
        if( base < 1 | base > n ){ stop( "baseline group number out of range" ) }
        contr <- contr[, , drop = FALSE]
    }
    contr
}
