\name{display}
\docType{genericFunction}
\alias{display}
\alias{display,lm-method}
\alias{display,bayesglm-method}
\alias{display,bayesglm.h-method}
\alias{display,glm-method}
\alias{display,mer-method}
\alias{display,polr-method}

\title{Functions for Processing lm, glm, mer and polr Output}

\description{This generic function gives a
    clean printout of lm, glm, mer, and polr objects.}

\usage{
display (object, ...)

\S4method{display}{lm}(object, digits=2, detail=FALSE)
\S4method{display}{bayesglm}(object, digits=2, detail=FALSE)
\S4method{display}{bayesglm.h}(object, digits=2, detail=FALSE)
\S4method{display}{glm}(object, digits=2, detail=FALSE)
\S4method{display}{mer}(object, digits=2, detail=FALSE)
\S4method{display}{polr}(object, digits=2, detail=FALSE)
}

\arguments{
  \item{object}{The output of a call to lm, glm, mer, polr, or related regressions 
    function with n data points and k predictors.}
  \item{...}{further arguments passed to or from other methods.}
  \item{digits}{number of significant digits to display.}
  \item{detail}{defaul is \code{FALSE}, if \code{TRUE}, display p-values or z-values}
  }

\details{This generic function gives a
    clean printout of lm, glm, mer and polr objects, focusing 
    on the most pertinent pieces of information: 
    the coefficients and their standard errors, the sample size, 
    number of predictors, residual standard deviation, and R-squared.
    Note: R-squared is automatically displayed to 2 digits, and deviances are
    automatically displayed to 1 digit, no matter what.
}

\value{Coefficients
    and their standard errors, the sample size, number of predictors, 
    residual standard
    deviation, and R-squared}
    
\references{Andrew Gelman and Jennifer Hill, 
            Data Analysis Using Regression and Multilevel/Hierarchical Models, 
            Cambridge University Press, 2006.}

\author{Andrew Gelman \email{gelman@stat.columbia.edu}; 
        Yu-Sung Su \email{ys463@columbia.edu};
        Maria Grazia Pittau \email{grazia@stat.columbia.edu}
        }
        
\note{Output are the model, the regression coefficients and standard errors,
 and the residual sd and R-squared (for a linear model),
 or the null deviance and residual deviance (for a generalized linear model).
}

\seealso{\code{\link[base]{summary}},
         \code{\link[stats]{lm}},
         \code{\link[stats]{glm}},
         \code{\link[lme4]{lmer}},
         \code{\link[MASS]{polr}}
         }

\examples{                          
# Here's a simple example of a model of the form, y = a + bx + error, 
# with 10 observations in each of 10 groups, and with both the 
# intercept and the slope varying by group.  First we set up the model and data.
   group <- rep(1:10, rep(10,10))
   group2 <- rep(1:10, 10)
   mu.a <- 0
   sigma.a <- 2
   mu.b <- 3
   sigma.b <- 4
   rho <- 0.56
   Sigma.ab <- array (c(sigma.a^2, rho*sigma.a*sigma.b, 
                    rho*sigma.a*sigma.b, sigma.b^2), c(2,2))
   sigma.y <- 1
   ab <- mvrnorm (10, c(mu.a,mu.b), Sigma.ab)
   a <- ab[,1]
   b <- ab[,2]
   d <- rnorm(10)

   x <- rnorm (100)
   y1 <- rnorm (100, a[group] + b*x, sigma.y)
   y2 <- rbinom(100, 1, prob=invlogit(a[group] + b*x))
   y3 <- rnorm (100, a[group] + b[group]*x + d[group2], sigma.y)
   y4 <- rbinom(100, 1, prob=invlogit(a[group] + b*x + d[group2]))


# display a simple linear model

   M1 <- lm (y1 ~ x)
   display (M1)
  

# display a simple logit model

   M2 <- glm (y2 ~ x, family=binomial(link="logit"))
   display (M2)

# Then fit and display a simple varying-intercept model:
 
   M3 <- lmer (y1 ~ x + (1|group))
   display (M3)
#   M3.sim <- mcsamp (M3)
#   print (M3.sim)
#   plot (M3.sim)
 
# Then the full varying-intercept, varying-slope model:
 
   M4 <- lmer (y1 ~ x + (1 + x |group))
   display (M4)
#   M4.sim <- mcsamp (M4)
#   print (M4.sim)
#   plot (M4.sim)
 
# Then the full varying-intercept, logit model:
 
   M5 <- glmer (y2 ~ x + (1|group), family=binomial(link="logit"))
   display (M5)
#   M5.sim <- mcsamp (M5)
#   print (M5.sim)
#   plot (M5.sim)
 
# Then the full varying-intercept, varying-slope logit model:
 
   M6 <- glmer (y2 ~ x + (1|group) + (0 + x |group), 
        family=binomial(link="logit"))
   display (M6)
#   M6.sim <- mcsamp (M6)
#   print (M6.sim)
#   plot (M6.sim)

# Then non-nested varying-intercept, varying-slop model:
   
   M7 <- lmer (y3 ~ x + (1 + x |group) + (1|group2))
   display(M7)
#   M7.sim <- mcsamp (M7)
#   print (M7.sim)
#   plot (M7.sim)

 
# Then the ordered logit model from polr

   M8 <- polr(Sat ~ Infl + Type + Cont, weights = Freq, data = housing)
   display(M8)
   
   M9 <- bayespolr(Sat ~ Infl + Type + Cont, weights = Freq, data = housing)
   display(M9) 
   

}


\keyword{manip}
\keyword{methods}
