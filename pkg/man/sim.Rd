\name{sim}
\docType{genericFunction}
\alias{sim}
\alias{sim,lm-method}
\alias{sim,glm-method}
\alias{sim,mer-method}
\alias{sim,lmer2-method}


\title{Functions to Get Posterior Distributions}

\description{
This generic function gets posterior simulations of sigma and beta from a lm object, or
simulations of beta from a glm object, or
simulations of beta from a mer object
}

\usage{
sim(object, ...)

\S4method{sim}{lm}(object, n.sims = 100)
\S4method{sim}{glm}(object, n.sims = 100)
\S4method{sim}{mer}(object, n.sims = 100)
\S4method{sim}{lmer2}(object, n.sims = 100)
}

\arguments{
  \item{object}{the output of a call to "lm" with n data points and k predictors.}
  \item{...}{further arguments passed to or from other methods.}
  \item{n.sims}{number of independent simulation draws to create.}

}

\value{
   \item{sigma.sim}{vector of n.sims random draws of sigma
  (for glm's, this just returns a vector of 1's or else of the
  square root of the overdispersion parameter if that is in the model)}
  \item{beta.sim}{matrix (dimensions n.sims x k) of n.sims random draws of beta}
  }

\references{Andrew Gelman and Jennifer Hill, 
            Data Analysis Using Regression and Multilevel/Hierarchical Models, 
            Cambridge University Press, 2006.}

\author{Andrew Gelman \email{gelman@stat.columbia.edu};
        Yu-Sung Su \email{ys463@columbia.edu};
        M.Grazia Pittau \email{grazia@stat.columbia.edu}
        }

\seealso{\code{\link{display}},
        \code{\link{mcsamp}},
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

# Example 2
 u.full <- u[group]
 E2 <- lmer (y ~ x + u.full + (1 | group))
 display(E2)
 E2.sim <- sim (E2)

# Example 3 
 y <- rbinom (n, 1, invlogit (a[group] + b*x))
 E3 <- lmer (y ~ x + (1 | group), family=binomial(link="logit"),
    control=list(usePQL=TRUE))
 display(E3)
 E3.sim <- sim (E3)
 
\testonly{
 E4 <- lmer2 (y ~ x + u.full + (1 | group))
 E4.sim <- sim (E4) 
}
}

\keyword{models}
\keyword{methods}
