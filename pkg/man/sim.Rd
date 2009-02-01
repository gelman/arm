\name{sim}
\docType{genericFunction}
\alias{sim}
\alias{sim,lm-method}
\alias{sim,glm-method}
\alias{sim,mer-method}



\title{Functions to Get Posterior Distributions}

\description{
This generic function gets posterior simulations of sigma and beta from a \code{lm} object, or
simulations of beta from a \code{glm} object, or
simulations of beta from a \code{mer} object
}

\usage{
sim(object, ...)

\S4method{sim}{lm}(object, n.sims = 100)
\S4method{sim}{glm}(object, n.sims = 100)
\S4method{sim}{mer}(object, n.sims = 100, ranef=TRUE)

}

\arguments{
  \item{object}{the output of a call to \code{lm} with n data points and k predictors.}
  \item{...}{further arguments passed to or from other methods.}
  \item{n.sims}{number of independent simulation draws to create.}
  \item{ranef}{default is \code{TRUE}, which will return simulated random effects.}
}

\value{
  \item{coef}{matrix (dimensions n.sims x k) of n.sims random draws of coefficients.
    Previously, it is called \code{beta}.}
  \item{fixef}{matrix (dimensions n.sims x k) of n.sims random draws of coefficients 
    of the fixed effects for the \code{mer} objects. Previously, it is called \code{unmodeled}.}
  \item{sigma}{vector of n.sims random draws of sigma
  (for \code{glm}'s, this just returns a vector of 1's or else of the
  square root of the overdispersion parameter if that is in the model)}
  }

\references{Andrew Gelman and Jennifer Hill. (2006).
            \emph{Data Analysis Using Regression and Multilevel/Hierarchical Models}.
            Cambridge University Press.}

\author{Andrew Gelman \email{gelman@stat.columbia.edu};
        Yu-Sung Su \email{ys463@columbia.edu}
}

\seealso{\code{\link{display}},
        \code{\link{lm}},
        \code{\link{glm}},
        \code{\link[lme4]{lmer}}
}

\examples{
#Examples of "sim" 
 set.seed (1)
 J <- 15
 n <- J*(J+1)/2
 group <- rep (1:J, 1:J)
 mu.a <- 5
 sigma.a <- 2
 a <- rnorm (J, mu.a, sigma.a)
 b <- -3
 x <- rnorm (n, 2, 1)
 sigma.y <- 6
 y <- rnorm (n, a[group] + b*x, sigma.y)
 u <- runif (J, 0, 3)
 y123.dat <- cbind (y, x, group)
# Linear regression 
 x1 <- y123.dat[,2]
 y1 <- y123.dat[,1]
 M1 <- lm (y1 ~ x1)
 display(M1)
 M1.sim <- sim(M1)
 ## to get the uncertainty for the simulated estimates
 apply(M1.sim$coef, 2, quantile)
 quantile(M1.sim$sigma)
 
# Logistic regression 
 u.data <- cbind (1:J, u)
 dimnames(u.data)[[2]] <- c("group", "u")
 u.dat <- as.data.frame (u.data)
 y <- rbinom (n, 1, invlogit (a[group] + b*x))
 M2 <- glm (y ~ x, family=binomial(link="logit"))
 display(M2)
 M2.sim <- sim (M2)

# Using lmer:
# Example 1
 E1 <- lmer (y ~ x + (1 | group))
 display(E1)
 E1.sim <- sim (E1)
 ## obtain intervals of the estimates
 t(apply(E1.sim$fixef, 2, quantile, probs = c(0.025, 0.5, 0.975)))


# Example 2
 u.full <- u[group]
 E2 <- lmer (y ~ x + u.full + (1 | group))
 display(E2)
 E2.sim <- sim (E2)

# Example 3 
 y <- rbinom (n, 1, invlogit (a[group] + b*x))
 E3 <- glmer (y ~ x + (1 | group), family=binomial(link="logit"))
 display(E3)
 E3.sim <- sim (E3)
 
}

\keyword{models}
\keyword{methods}
