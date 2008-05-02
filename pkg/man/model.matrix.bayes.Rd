\name{model.matrix.bayes}
%\docType{genericFunction}
\alias{model.matrix.bayes}
\alias{model.matrix.bayes.h}
%\alias{model.matrix.bayes,bayesglm-method}
%\alias{model.matrix.bayes,bayesglm.h-method}

\title{Construct Design Matrices}

\description{
  \code{model.matrix.bayes} creates a design matrix.
}

\usage{
model.matrix.bayes(object, data = environment(object), 
    contrasts.arg = NULL, xlev = NULL, keep.order = FALSE, drop.baseline=FALSE,...)
    
model.matrix.bayes.h(object, data = environment(object), 
    contrasts.arg = NULL, xlev = NULL, keep.order = FALSE,  batch = NULL, ...)
}

\arguments{
  \item{object}{an object of an appropriate class.  For the default
    method, a model formula or terms object.}
  \item{data}{a data frame created with \code{\link{model.frame}}.  If
    another sort of object, \code{model.frame} is called first.}
  \item{contrasts.arg}{A list, whose entries are contrasts suitable for
    input to the \code{\link{contrasts}} replacement function and
    whose names are the names of columns of \code{data} containing
    \code{\link{factor}}s.}
  \item{xlev}{to be used as argument of \code{\link{model.frame}} if
    \code{data} has no \code{"terms"} attribute.}
  \item{keep.order}{a logical value indicating whether the terms should
    keep their positions. If \code{FALSE} the terms are reordered so
    that main effects come first, followed by the interactions,
    all second-order, all third-order and so on.  Effects of a given
    order are kept in the order specified.} 
  \item{drop.baseline}{Drop the base level of categorical x's, default is TRUE.}
  \item{batch}{Not implement yet!}
  \item{\dots}{further arguments passed to or from other methods.}
}
\details{
\code{model.matrix.bayes} is adapted from \code{model.matrix} in the \code{stats}
pacakge and is designed for the use of \code{\link{bayesglm}} and \code{bayesglm.hierachical} (not yet implemented!).
It is designed to keep baseline levels of all categorical varaibles and keep the
variable names unodered in the output.  The design matrices created by 
\code{model.matrix.bayes} are unidentifiable using classical regression methods, 
though; they can be identified using \code{\link{bayesglm}} and 
\code{bayesglm.hierachical}.
}

\references{Andrew Gelman, Aleks Jakulin, Maria Grazia Pittau and Yu-Sung Su, 
            A default prior distribution for logistic and other regression models,  
unpublished paper available at \url{http://www.stat.columbia.edu/~gelman/standardize/}
}

\seealso{
  \code{\link{model.frame}}, \code{\link{model.extract}},
  \code{\link{terms}}, \code{\link{terms.formula}}, \code{\link{bayesglm}}. 
}

\author{Yu-Sung Su \email{ys463@columbia.edu}}

\examples{
ff <- log(Volume) ~ log(Height) + log(Girth)
str(m <- model.frame(ff, trees))
(model.matrix(ff, m))
class(ff) <- c("bayesglm", "terms", "formula")
(model.matrix.bayes(ff, m))
class(ff) <- c("bayesglm.h", "terms", "formula")
(model.matrix.bayes(ff, m))
}
\keyword{models}
