\name{sigma.hat}
%\docType{genericFunction}
\alias{sigma.hat}
\alias{sigma.hat,lm-method}
\alias{sigma.hat,glm-method}
\alias{sigma.hat,merMod-method}
\alias{sigma.hat,sim-method}
\alias{sigma.hat,sim.merMod-method}



\title{Extract Residual Errors}

\description{This generic function extracts residual errors from a fitted model.
}

\usage{
sigma.hat(object,\dots)

\S4method{sigma.hat}{lm}(object)
\S4method{sigma.hat}{glm}(object)
\S4method{sigma.hat}{merMod}(object)
\S4method{sigma.hat}{sim}(object)
\S4method{sigma.hat}{sim.merMod}(object)

}

\arguments{
  \item{object}{any fitted model object of \code{lm}, \code{glm} and \code{mer} class}
  \item{\dots}{other arguments}
}

\author{Andrew Gelman \email{gelman@stat.columbia.edu};
        Yu-Sung Su \email{suyusung@tsinghua.edu.cn}
    }

\seealso{\code{\link{display}},
        \code{\link{summary}},
        \code{\link{lm}},
        \code{\link{glm}},
        \code{\link[lme4]{lmer}}
}

\examples{
   group <- rep(1:10, rep(10,10))
   mu.a <- 0
   sigma.a <- 2
   mu.b <- 3
   sigma.b <- 4
   rho <- 0
   Sigma.ab <- array (c(sigma.a^2, rho*sigma.a*sigma.b, 
                    rho*sigma.a*sigma.b, sigma.b^2), c(2,2))
   sigma.y <- 1
   ab <- mvrnorm (10, c(mu.a,mu.b), Sigma.ab)
   a <- ab[,1]
   b <- ab[,2]

   x <- rnorm (100)
   y1 <- rnorm (100, a[group] + b[group]*x, sigma.y)
   y2 <- rbinom(100, 1, prob=invlogit(a[group] + b*x))
   
   M1 <- lm (y1 ~ x) 
   sigma.hat(M1)
   
   M2 <- bayesglm (y1 ~ x, prior.scale=Inf, prior.df=Inf)
   sigma.hat(M2) # should be same to sigma.hat(M1)

   M3 <- glm (y2 ~ x, family=binomial(link="logit"))
   sigma.hat(M3)

   M4 <- lmer (y1 ~ (1+x|group))
   sigma.hat(M4) 
   
   M5 <- glmer (y2 ~ (1+x|group), family=binomial(link="logit"))
   sigma.hat(M5)

}
\keyword{manip}
\keyword{methods}
