\name{extractDIC}
%\docType{genericFunction}
\alias{extractDIC}
\alias{extractDIC.merMod}
%\alias{extractAIC.mer}

\title{
Extract DIC from a \code{merMod} model
}

\description{
 Computes the (generalized) *D*eviance *I*nformation *C*riterion for a \code{merMod} model.
}
\usage{
extractDIC(fit,\dots)
\S3method{extractDIC}{merMod}(fit,\dots)

%\method{extractAIC}{mer}(fit,\dots)
}

\arguments{
  \item{fit}{fitted \code{merMod} mode, usually the result of a fiiter like \code{merMod}.}
  \item{\dots}{further arguments (currently unused).}
}


\author{
  Andrew Gelman \email{gelman@stat.columbia.edu};
  Yu-Sung Su \email{suyusung@tsinghua.edu.cn}
}


\examples{
fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
AIC(fm1)
extractDIC(fm1)
}

\keyword{manip}
\keyword{methods}
