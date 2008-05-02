# ========================================================================
# function for trace plot
# ========================================================================
traceplot.bugs <- function( x, mfrow = c( 1, 1 ), varname = NULL,
  match.head = TRUE, ask = TRUE,
  col = rainbow( x$n.chains ),
  lty = 1, lwd = 1, ... ) 
{
  par( mfrow = mfrow )
  par( ask = ask )
  n.chain    <- x$n.chains
  n.keep     <- x$n.keep
  bugs.array <- x$sims.array
  varnamelist <- gsub( "\\[.*\\]","", dimnames( bugs.array )[[3]], fixed = FALSE )
  if( is.null( varname ) ){ varname <- ".*" }
  if( match.head ) { varname <- paste( "^", varname, sep="" ) }
  index      <- unlist( sapply( varname, function( x ){ grep( x, varnamelist ) } ) )
  n.var      <- length( index )
  for( j in index ) {
    range.x  <- c( 1, n.keep )
    range.y  <- range( bugs.array[,,j] )
    v.name   <- dimnames( bugs.array )[[3]][j]
    plot( range.x, range.y, type = "n", main = v.name,
            xlab = "iteration", ylab = v.name,
              xaxt  = "n", xaxs = "i", ... )
    for( i in 1:n.chain ) {
      x.cord <- 1:n.keep
      y.cord <- bugs.array[,i,j]
      lines( x.cord , y.cord , col = col[i], lty = lty, lwd = lwd )
    }
    axis( 1, at = seq(0, n.keep, n.keep*0.1), tick = TRUE )
  }
}
