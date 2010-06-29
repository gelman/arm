\name{extractDIC}
\docType{genericFunction}
\alias{extractDIC}
\alias{extractDIC.mer}
\alias{extractAIC.mer}

\title{
Extract AIC and DIC from a mer model
}

\description{
 Computes the (generalized) Akaike *A*n *I*nformation *C*riterion 
 and *D*eviance *I*nformation *C*riterion for a mer model.
}
\usage{
extractDIC(fit,\dots)
\method{extractDIC}{mer}(fit,\dots)

\method{extractAIC}{mer}(fit,\dots)
}

\arguments{
  \item{fit}{fitted \code{mer} mode, usually the result of a fiiter like \code{mer}.}
  \item{\dots}{further arguments (currently unused).}
}


\author{
  Andrew Gelman \email{gelman@stat.columbia.edu};
  Yu-Sung Su \email{ys463@columbia.edu}
}


\examples{
fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
extractAIC(fm1)
extractDIC(fm1)
}

\keyword{manip}
\keyword{methods}
