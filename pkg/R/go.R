#' Function to Recall Last Source File
#'
#' A function that like \code{source()} but recalls the last source
#' file names by default.
#' 
#' @param \dots list of filenames as character strings.
#' @param add add these names to the current list; if replace, then
#' \code{FALSE}.
#' @param timer time the execution time of \code{go()}.
#'
#' @author Jouni Kerman \email{jouni@@kerman.com}
#' \email{kerman@@stat.columbia.edu}
#' @keywords methods manip
#' @export
#' @examples
#' go('myprog')            # will run source('myprog.r')
#' go()                    # will run source('myprog.r') again
#' go('somelib',add=TRUE)  # will run source('myprog.r') and source('somelib.r')
#' go('myprog','somelib')  # same as above
#' go('mytest')            # will run source('mytest') only
#' go()                    # runs source('mytest') again
#' G                       # short cut to call go()
go <- function(..., add=FALSE, timer=FALSE)
# Modified:    2004-06-22
{
  last.sources <- getOption(".Last.Source")
  sources <- unlist(list(...))
  if (length(sources)<1) {
    sources <- last.sources
  } else if (add) {
    sources <- c(last.sources,sources)
  }
  if (length(sources)<1) {
    return(cat("Usage: go('sourcefile', 'sourcefile2', ..., add=?, timer=?)\n"))
  }
  options(".Last.Source"=sources)
  cat("Source file(s): ",sources,"\n")
  yy <- NULL
  for (src in sources) {
    if (is.na(src)) {
      next
    }
    if (!file.exists(src)) {
      src2 <- paste(src, ".R", sep="")
      if (file.exists(src2))
        src <- src2
      else {
        cat("source('",src,"') : file does not exist.\n",sep='')
        next
      }
    }
    cat("source('",src,"')\n",sep="")
    if (timer)
      cat("source('",src,"') : ",max(na.omit(system.time(source(src)))),
        " seconds elapsed.\n", sep='')
    else 
      yy[[src]] <- source(src)
  }
  invisible(yy)
}


# By entering "G" on the console, go() is run. This is faster than
# typing "go()"...
#' @rdname go
#' @method print GO
#' @export
print.GO <- function(x,...) {go()}

#' @rdname go
#' @export
G <- structure(NA, class="GO")
#class(G) <- "GO"

# end of go.R
