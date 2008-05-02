\name{standardize}
\alias{standardize}

\title{Function for Standardizing Regression Predictors by Centering and 
    Dividing by 2 sd's}

\description{Numeric variables that take on more than two values are each rescaled
     to have a mean of 0 and a sd of 0.5;
     Binary variables are rescaled to have a mean of 0 and a difference of 1
     between their two categories; 
     Non-numeric variables that take on more than two values are unchanged;
     Variables that take on only one value are unchanged

}
\usage{
standardize(object, unchanged = NULL, 
    standardize.y = FALSE, binary.inputs = "center")
}

\arguments{
  \item{object}{an object of class "lm" or "glm"}
  \item{unchanged}{vector of names of parameters to leave unstandardized}
  \item{standardize.y}{ if TRUE, the outcome variable is standardized also}
  \item{binary.inputs}{options for standardizing binary variables}
}

\details{
    "0/1" (rescale so that the lower value is 0 and the upper is 1)
    "-0.5/0.5" (rescale so that the lower value is -0.5 and upper is 0.5)
    "center" (rescale so that the mean of the data is 0 and the difference
              between the two categories is 1)
    "full" (rescale by subtracting the mean and dividing by 2 sd's)
    "leave.alone" (do nothing)
}

\references{Andrew Gelman, 
    Scaling regression inputs by dividing by two standard deviations 
    \url{http://www.stat.columbia.edu/~gelman/research/unpublished/standardizing.pdf}
}

\author{Andrew Gelman \email{gelman@stat.columbia.edu}
        Yu-Sung Su \email{ys463@columbia.edu}
}

\seealso{\code{\link{rescale}}
            }



\examples{
  # Set up the fake data
  n <- 100
  x <- rnorm (n, 2, 1)
  x1 <- rnorm (n)
  x1 <- (x1-mean(x1))/(2*sd(x1))   # standardization
  x2 <- rbinom (n, 1, .5)
  b0 <- 1
  b1 <- 1.5
  b2 <- 2
  y <- rbinom (n, 1, invlogit(b0+b1*x1+b2*x2))
  M1 <- glm (y ~ x, family=binomial(link="logit"))
  display(M1)
  M2 <-standardize(M1)
  display(M2)
  }
\keyword{models}
\keyword{methods}
