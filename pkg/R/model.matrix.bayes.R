#' Construct Design Matrices
#'
#' \code{model.matrixBayes} creates a design matrix.
#'
#' \code{model.matrixBayes} is adapted from \code{model.matrix} in the
#' \code{stats} pacakge and is designed for the use of
#' \code{\link{bayesglm}}.% and \code{bayesglm.hierachical} (not yet
#' implemented!). It is designed to keep baseline levels of all
#' categorical varaibles and keep the variable names unodered in the
#' output.  The design matrices created by \code{model.matrixBayes}
#' are unidentifiable using classical regression methods, though; they
#' can be identified using \code{\link{bayesglm}}.% and
#' %\code{bayesglm.hierachical}.
#' 
#' @param object an object of an appropriate class.  For the default
#'   method, a model formula or terms object.
#' @param data a data frame created with \code{\link{model.frame}}.  If
#'   another sort of object, \code{model.frame} is called first.
#' @param contrasts.arg A list, whose entries are contrasts suitable for
#'   input to the \code{\link{contrasts}} replacement function and
#'   whose names are the names of columns of \code{data} containing
#'   \code{\link{factor}}s.
#' @param xlev to be used as argument of \code{\link{model.frame}} if
#'   \code{data} has no \code{"terms"} attribute.
#' @param keep.order a logical value indicating whether the terms should
#'   keep their positions. If \code{FALSE} the terms are reordered so
#'   that main effects come first, followed by the interactions,
#'   all second-order, all third-order and so on.  Effects of a given
#'   order are kept in the order specified.
#' @param drop.baseline Drop the base level of categorical Xs,
#'   default is TRUE.
#' @param ... further arguments passed to or from other methods.
#' @references Andrew Gelman, Aleks Jakulin, Maria Grazia Pittau and
#' Yu-Sung Su. (2009). \dQuote{A Weakly Informative Default Prior
#' Distribution For Logistic And Other Regression Models.} \emph{The
#' Annals of Applied Statistics} 2 (4):
#' 1360--1383. \url{http://www.stat.columbia.edu/~gelman/research/published/priors11.pdf}
#' @author Yu-Sung Su \email{suyusung@@tsinghua.edu.cn}
#' @seealso \code{\link[stats]{model.frame}},
#' \code{\link[stats]{model.extract}}, \code{\link[stats]{terms}},
#' \code{\link[stats]{terms.formula}}, \code{\link{bayesglm}}.
#' @keywords models manip
#' @export
#' @examples
#' ff <- log(Volume) ~ log(Height) + log(Girth)
#' str(m <- model.frame(ff, trees))
#' (model.matrix(ff, m))
#' class(ff) <- c("bayesglm", "terms", "formula")
#' (model.matrixBayes(ff, m))
model.matrixBayes <- function(object, data = environment(object),
        contrasts.arg = NULL, xlev = NULL, keep.order=FALSE, drop.baseline=FALSE,...)
{
    #class(object) <- c("terms", "formula")
    t <- if( missing( data ) ) { 
          terms( object ) 
         }else{ 
            terms.formula(object, data = data, keep.order=keep.order) 
         }
    attr(t, "intercept") <- attr(object, "intercept")
    if (is.null(attr(data, "terms"))){ 
      data <- model.frame(object, data, xlev=xlev) 
    }else {
        reorder <- match(sapply(attr(t,"variables"), deparse, width.cutoff=500)[-1], names(data))
        if (any(is.na(reorder))) {
            stop( "model frame and formula mismatch in model.matrix()" ) 
        }
        if(!identical(reorder, seq_len(ncol(data)))) {
            data <- data[,reorder, drop = FALSE] 
        }
    }
    int <- attr(t, "response")
    if(length(data)) {      # otherwise no rhs terms, so skip all this
        
        if (drop.baseline){
          contr.funs <- as.character(getOption("contrasts"))
        }else{
          contr.funs <- as.character(list("contr.bayes.unordered", "contr.bayes.ordered"))
        }
        
        namD <- names(data)
        ## turn any character columns into factors
        for(i in namD)
            if(is.character( data[[i]] ) ) {
                data[[i]] <- factor(data[[i]])
                warning( gettextf( "variable '%s' converted to a factor", i ), domain = NA)
            }
        isF <- vapply(data, function(x) is.factor(x) || is.logical(x), NA)        
        isF[int] <- FALSE
        isOF <- vapply(data, is.ordered, NA)
        for( nn in namD[isF] )            # drop response
            if( is.null( attr( data[[nn]], "contrasts" ) ) ) {
                contrasts( data[[nn]] ) <- contr.funs[1 + isOF[nn]]
            }
        ## it might be safer to have numerical contrasts:
        ##    get(contr.funs[1 + isOF[nn]])(nlevels(data[[nn]]))
        if ( !is.null( contrasts.arg ) && is.list( contrasts.arg ) ) {
            if ( is.null( namC <- names( contrasts.arg ) ) ) {
                stop( "invalid 'contrasts.arg' argument" )
            }
            for (nn in namC) {
                if ( is.na( ni <- match( nn, namD ) ) ) {
                    warning( gettextf( "variable '%s' is absent, its contrast will be ignored", nn ), domain = NA )
                }
                else {
                    ca <- contrasts.arg[[nn]]
                    if( is.matrix( ca ) ) {
                        contrasts( data[[ni]], ncol( ca ) ) <- ca
                    }
                    else { 
                        contrasts( data[[ni]] ) <- contrasts.arg[[nn]]
                    }
                }
            }
        }
    } else {               # internal model.matrix needs some variable
        isF  <-  FALSE
        data <- data.frame(x=rep(0, nrow(data)))
    }
    #ans  <- .Internal( model.matrix( t, data ) )
    ans  <- model.matrix.default(object=t, data=data)
    cons <- if(any(isF)){
              lapply( data[isF], function(x) attr( x,  "contrasts") ) 
            }else { NULL }
    attr(ans, "contrasts" ) <- cons
    ans
}


#setMethod("model.matrix.bayes", signature(object = "bayesglm.h"),
#model.matrix.bayes.h <- function (object, data = environment(object), 
#            contrasts.arg = NULL, 
#            xlev = NULL, keep.order = FALSE, batch = NULL, ...)
#{
#    class(object) <- c("formula")
#    t <- if (missing(data)) {
#        terms(object)
#    }
#    else {
#        terms(object, data = data, keep.order = keep.order)
#    }
#    attr(t, "intercept") <- attr(object, "intercept")
#    if (is.null(attr(data, "terms"))) {
#        data <- model.frame(object, data, xlev = xlev)
#    }
#    else {
#        reorder <- match(sapply(attr(t, "variables"), deparse, 
#            width.cutoff = 500)[-1], names(data))
#        if (any(is.na(reorder))) {
#            stop("model frame and formula mismatch in model.matrix()")
#        }
#        if (!identical(reorder, seq_len(ncol(data)))) {
#            data <- data[, reorder, drop = FALSE]
#        }
#    }
#    int <- attr(t, "response")
#    if (length(data)) {
#        contr.funs <- as.character(getOption("contrasts"))
#        contr.bayes.funs <- as.character(list("contr.bayes.unordered", 
#            "contr.bayes.ordered"))
#        namD <- names(data)
#        for (i in namD) if (is.character(data[[i]])) {
#            data[[i]] <- factor(data[[i]])
#            warning(gettextf("variable '%s' converted to a factor", i), domain = NA)
#        }
#        isF <- sapply(data, function(x) is.factor(x) || is.logical(x))
#        isF[int] <- FALSE
#        isOF <- sapply(data, is.ordered)
#        if (length(batch) > 1) {
#            ba <- batch[isF[-1]]
#        }
#        else if (length(batch) == 1) {
#            ba <- rep(batch, length(isF[-1]))
#        }
#        else {
#            ba <- rep(0, length(isF[-1]))
#        }
#        iin <- 1
#        for (nn in namD[isF]) if (is.null(attr(data[[nn]], "contrasts"))) {
#            if (ba[[iin]] > 0) {
#                contrasts(data[[nn]]) <- contr.bayes.funs
#            }
#            else {
#                contrasts(data[[nn]]) <- contr.funs
#            }
#            iin <- iin + 1
#        }
#        if (!is.null(contrasts.arg) && is.list(contrasts.arg)) {
#            if (is.null(namC <- names(contrasts.arg))) {
#                stop("invalid 'contrasts.arg' argument")
#            }
#            for (nn in namC) {
#                if (is.na(ni <- match(nn, namD))) {
#                  warning(gettextf("variable '%s' is absent, its contrast will be ignored", 
#                    nn), domain = NA)
#                }
#                else {
#                  ca <- contrasts.arg[[nn]]
#                  if (is.matrix(ca)) {
#                    contrasts(data[[ni]], ncol(ca)) <- ca
#                  }
#                  else {
#                    contrasts(data[[ni]]) <- contrasts.arg[[nn]]
#                  }
#                }
#            }
#        }
#    }
#    else {
#        isF <- FALSE
#        data <- list(x = rep(0, nrow(data)))
#    }
#    ans <- .Internal(model.matrix(t, data))
#    cons <- if (any(isF)) {
#        lapply(data[isF], function(x) attr(x, "contrasts"))
#    }
#    else {
#        NULL
#    }
#    attr(ans, "contrasts") <- cons
#    ans
#}
##)
