#' residual plot for the observed values
#'
#' Plots the residual of observed variable.
#'
#' @param Expected Expected value.
#' @param Residuals Residual value.
#' @param sigma Standard error.
#' @param main main for the plot. See \code{plot} for detail.
#' @param col.pts Color of the points.
#' @param col.ctr Color of the line at zero.
#' @param col.sgm Color of standard error line.
#' @param cex A numerical value giving the amount by which plotting
#' text and symbols should be magnified relative to the default.  See
#' par for detail.
#' @param gray.scale If \code{TRUE}, makes the plot into black and
#' white. This option overwrites the color specification. Default is
#' FALSE.
#' @param xlab Label for x axis
#' @param ylab Label for y axis
#' @param ... Additional parameters passed to \code{plot} function
#' @return Plot to visualize pattern of residulal value for the
#' expected value.
#' @author 
#'   Masanao Yajima \email{yajima@@stat.columbia.edu},
#'   M.Grazia Pittau \email{grazia@@stat.columbia.edu}
#' @keywords hplot
#' @export
#' @examples
#' old.par <- par(no.readonly = TRUE)
#' 
#' x <- rnorm(100)
#' y <- rnorm(100)
#' fit <- lm(y~x)
#' y.hat <- fitted(fit)
#' u <- resid(fit)
#' sigma <- sigma.hat(fit)
#' residual.plot(y.hat, u, sigma)
#' 
#' par(old.par)




residual.plot <- function ( Expected, Residuals, sigma, 
                            main = deparse(substitute( Expected )), 
                            col.pts = "blue", col.ctr = "red", 
                            col.sgm = "black", cex = 0.5, gray.scale = FALSE, 
                            xlab="Predicted", ylab="Residuals", ... ) {
  if( gray.scale == TRUE ) { 
    col.pts <- "black";
    col.ctr <- "black";
    col.sgm <- "gray60";
  }
  plot( Expected[!is.na( Residuals )], Residuals[ !is.na( Residuals ) ],
         xlab = xlab, ylab = ylab, main = main, col = col.pts,
          pch = 19, cex = cex, ... );
  #mtext( "Residuals vs Predicted", 3, cex= 0.6 )  #, adj=1 );
  # add the zero line for clarity
  abline ( h = 0, lty = "dashed", col = col.ctr );
  # residual s.e.
  resid.se <- sigma;
  # Add two-standard-error lines
  abline ( h =  2*resid.se, lty = "dashed", col = col.sgm );
  abline ( h = -2*resid.se, lty = "dashed", col = col.sgm );
} 
