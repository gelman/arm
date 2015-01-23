#' Generic Function for Making Coefficient Plot
#'
#' Functions that plot the coefficients plus and minus 1 and 2 sd from
#' a lm, glm, bugs, and polr fits.
#'
#' This function plots coefficients from bugs, lm, glm and polr with 1
#' sd and 2 sd interval bars.
#'
#' @param object fitted objects-lm, glm, bugs and polr, or a vector of
#'   coefficients.
#' @param ... further arguments passed to or from other methods.
#' @param coefs a vector of coefficients.
#' @param sds a vector of sds of coefficients.
#' @param CI confidence interval, default is 2, which will plot plus
#'   and minus 2 sds or 95\% CI. If CI=1, plot plus and minus 1 sds or
#'   50\% CI instead.
#' @param lower.conf.bounds lower bounds of confidence intervals.
#' @param upper.conf.bounds upper bounds of confidence intervals.
#' @param varnames a vector of variable names, default is NULL, which
#'   will use the names of variables; if specified, the length of
#'   varnames must be equal to the length of predictors, including the
#'   intercept.
#' @param vertical orientation of the plot, default is TRUE which will
#'   plot variable names in the 2nd axis. If FALSE, plot variable names
#'   in the first axis instead.
#' @param v.axis default is TRUE, which shows the bottom
#'   axis--axis(1).
#' @param h.axis default is TRUE, which shows the left axis--axis(2).
#' @param cex.var The fontsize of the varible names, default=0.8.
#' @param cex.pts The size of data points, default=0.9.
#' @param col.pts color of points and segments, default is black.
#' @param pch.pts symbol of points, default is solid dot.
#' @param var.las the orientation of variable names against the axis,
#'   default is 2. See the usage of \code{las} in \code{\link{par}}.
#' @param main The main title (on top) using font and size (character
#'   expansion) \code{par("font.main")} and color
#'   \code{par("col.main")}.
#' @param xlab X axis label using font and character expansion
#'   \code{par("font.lab")} and color \code{par("col.lab")}.
#' @param ylab Y axis label, same font attributes as \code{xlab}.
#' @param mar A numerical vector of the form \code{c(bottom, left,
#'   top, right)} which gives the number of lines of margin to be
#'   specified on the four sides of the plot. The default is
#'   \code{c(1,3,5.1,2)}.
#' @param plot default is TRUE, plot the estimates.
#' @param add if add=TRUE, plot over the existing plot. default is
#'   FALSE.
#' @param offset add extra spaces to separate from the existing
#'   dots. default is 0.1.
## @param lower.bound default is -Inf.
#' @param var.idx the index of the variables of a bugs object, default
#'   is NULL which will plot all the variables.
#' @param intercept If TRUE will plot intercept, default=FALSE to get
#'   better presentation.
#'
#' @return Plot of the coefficients from a bugs, lm or glm fit. You
#'   can add the intercept, the variable names and the display the
#'   result of the fitted model.
#' @references Andrew Gelman and Jennifer Hill, Data Analysis Using
#'   Regression and Multilevel/Hierarchical Models, Cambridge University
#'   Press, 2006.
#' @author Yu-Sung Su \email{suyusung@@tsinghua.edu.cn}
#' @seealso \code{\link{display}}, \code{\link[graphics]{par}},
#'   \code{\link[stats]{lm}}, \code{\link[stats]{glm}},
#'   \code{\link{bayesglm}}, \code{\link[graphics]{plot}}
#' @keywords hplot dplot methods manip
#' @name coefplot
#' @examples
#' old.par <- par(no.readonly = TRUE)
#' 
#' y1 <- rnorm(1000,50,23)
#' y2 <- rbinom(1000,1,prob=0.72)
#' x1 <- rnorm(1000,50,2) 
#' x2 <- rbinom(1000,1,prob=0.63) 
#' x3 <- rpois(1000, 2) 
#' x4 <- runif(1000,40,100) 
#' x5 <- rbeta(1000,2,2) 
#' 
#' longnames <- c("a long name01","a long name02","a long name03",
#'                "a long name04","a long name05")
#'
#' fit1 <- lm(y1 ~ x1 + x2 + x3 + x4 + x5)
#' fit2 <- glm(y2 ~ x1 + x2 + x3 + x4 + x5, 
#'            family=binomial(link="logit"))
#' op <- par()
#' # plot 1
#' par (mfrow=c(2,2))
#' coefplot(fit1)
#' coefplot(fit2, col.pts="blue")
#' 
#' # plot 2
#' longnames <- c("(Intercept)", longnames) 
#' coefplot(fit1, longnames, intercept=TRUE, CI=1)
#' 
#' # plot 3
#' coefplot(fit2, vertical=FALSE, var.las=1, frame.plot=TRUE)
#' 
#' # plot 4: comparison to show bayesglm works better than glm
#' n <- 100
#' x1 <- rnorm (n)
#' x2 <- rbinom (n, 1, .5)
#' b0 <- 1
#' b1 <- 1.5
#' b2 <- 2
#' y <- rbinom (n, 1, invlogit(b0+b1*x1+b2*x2))
#' y <- ifelse (x2==1, 1, y)
#' x1 <- rescale(x1)
#' x2 <- rescale(x2, "center")
#' 
#' M1 <- glm (y ~ x1 + x2, family=binomial(link="logit"))
#'       display (M1)
#' M2 <- bayesglm (y ~ x1 + x2, family=binomial(link="logit"))
#'       display (M2)
#'
#' #=================== 
#' #    stacked plot
#' #===================
#' coefplot(M2, xlim=c(-1,5), intercept=TRUE)
#' coefplot(M1, add=TRUE, col.pts="red")
#'   
#' #==================== 
#' # arrayed plot       
#' #====================
#' par(mfrow=c(1,2))
#' x.scale <- c(0, 7.5) # fix x.scale for comparison
#' coefplot(M1, xlim=x.scale, main="glm", intercept=TRUE)
#' coefplot(M2, xlim=x.scale, main="bayesglm", intercept=TRUE)
#' 
#' # plot 5: the ordered logit model from polr
#' library("MASS")
#' M3 <- polr(Sat ~ Infl + Type + Cont, weights = Freq, data = housing)
#' coefplot(M3, main="polr")
#'   
#' M4 <- bayespolr(Sat ~ Infl + Type + Cont, weights = Freq, data = housing)
#' coefplot(M4, main="bayespolr", add=TRUE, col.pts="red")
#'
#' ## plot 6: plot bugs & lmer
#' # par <- op
#' # M5 <- lmer(Reaction ~ Days + (1|Subject), sleepstudy)
#' # M5.sim <- mcsamp(M5)
#' # coefplot(M5.sim, var.idx=5:22, CI=1, ylim=c(18,1), main="lmer model")
#'
#'
#' # plot 7: plot coefficients & sds vectors
#' coef.vect <- c(0.2, 1.4, 2.3, 0.5)
#' sd.vect <- c(0.12, 0.24, 0.23, 0.15)
#' longnames <- c("var1", "var2", "var3", "var4")
#' coefplot (coef.vect, sd.vect, varnames=longnames, main="Regression Estimates")
#' coefplot (coef.vect, sd.vect, varnames=longnames, vertical=FALSE, 
#'    var.las=1, main="Regression Estimates")
#'
#' par(old.par)
NULL

#' @rdname coefplot
#' @export
coefplot.default <- function(coefs, sds, CI=2, 
            lower.conf.bounds, upper.conf.bounds,
            varnames=NULL, 
            vertical=TRUE,
            v.axis=TRUE, h.axis=TRUE,
            cex.var=0.8, cex.pts=0.9, col.pts=1, pch.pts=20,
            var.las=2, main=NULL, xlab=NULL, ylab=NULL, mar=c(1,3,5.1,2),
            plot=TRUE, add=FALSE, offset=0.1, ...)
{
  
     # collect informations
    if (is.list(coefs)){
      coefs <- unlist(coefs)
    }
    
    
    n.x <- length(coefs)
    idx <- seq(1, n.x)   
    #bound <- lower.bound
    if(!missing(lower.conf.bounds)){
      if(length(coefs)!=length(lower.conf.bounds)){
        stop("Number of conf.bounds does not equal to number of estimates")
      }
    }
    if(!missing(upper.conf.bounds)){
      if(length(coefs)!=length(upper.conf.bounds)){
        stop("Number of conf.bounds does not equal to number of estimates")
      }
    }

    if(!missing(sds)){
      coefs.h <- coefs + CI*sds 
      coefs.l <- coefs - CI*sds
      est1 <- cbind(coefs - sds, coefs + sds)
      est2 <- cbind(coefs - 2*sds, coefs + 2*sds)
      if(!missing(lower.conf.bounds)){
        est1[,1] <- lower.conf.bounds
        CI <- 1
      }
      if(!missing(upper.conf.bounds)){
        est1[,2] <- upper.conf.bounds
        CI <- 1
      }

    }else{
      #coefs.h <- upper.conf.bounds
      #coefs.l <- lower.conf.bounds
      est1 <- cbind(coefs, coefs)
      if(!missing(lower.conf.bounds)){
        est1[,1] <- lower.conf.bounds
        CI <- 1
      }
      if(!missing(upper.conf.bounds)){
        est1[,2] <- upper.conf.bounds
        CI <- 1
      }
    }
    old.par <- par(no.readonly=TRUE)
    #on.exit(par(old.par))  
    min.mar <- par('mar')
    
    if (is.null(main)){main <- "Regression Estimates"}
    if (is.null(xlab)){xlab <- ""}
    if (is.null(ylab)){ylab <- ""}
        
    par(mar = mar)
    
    if (is.null(varnames)) {
      maxchar <- 0
    }
    else{
      maxchar <- max(sapply(varnames, nchar))
    }

    
    # add margin to the axis
    k <- 1/n.x   
    if(plot){
      if (vertical){

        mar[2] <- max(min.mar[2], trunc(mar[2] + maxchar/10)) + 0.1

        par(mar=mar)
        if(!add){
          plot(c(coefs.l, coefs.h), c(idx+k,idx-k), type="n",                                     
            axes=F, main=main, xlab=xlab, ylab=ylab,...) 
          if (h.axis){                                                  
            #axis(1)                                
            axis(3)
          }
          if (v.axis){
            axis(2, n.x:1, varnames[n.x:1], las=var.las, tck=FALSE, 
              lty=0, cex.axis=cex.var) 
          }
          abline(v=0, lty=2)                                                 
          points(coefs, idx, pch=pch.pts, cex=cex.pts, col=col.pts)
          if (CI==2){
            segments (est1[,1], idx, est1[,2], idx, lwd=2, col=col.pts)     
            segments (est2[,1], idx, est2[,2], idx, lwd=1, col=col.pts)
          }
          else{
            segments (est1[,1], idx, est1[,2], idx, lwd=1, col=col.pts)     
          }
        }
        else{
          idx <- idx + offset
          points(coefs, idx, pch=pch.pts, cex=cex.pts, col=col.pts)
          if (CI==2){
            segments (est1[,1], idx, est1[,2], idx, lwd=2, col=col.pts)     
            segments (est2[,1], idx, est2[,2], idx, lwd=1, col=col.pts)
          }
          else{
            segments (est1[,1], idx, est1[,2], idx, lwd=1, col=col.pts)     
          }
        }
    } # end of if vertical
    else{ # horizontal
      mar[1] <- max(min.mar[1], trunc(mar[1] + maxchar/10)) + 0.1
      par(mar=mar)
      if(!add){
        plot(c(idx+k,idx-k), c(coefs.l, coefs.h), type="n", axes=F, 
          main=main, xlab=xlab, ylab=ylab,...)                                                  
        if (v.axis){
          axis(2, las=var.las)                                
          #axis(4, las=var.las)
        }
        if (h.axis){
          axis(1, 1:n.x, varnames[1:n.x], las=var.las, tck=FALSE, 
            lty=0, cex.axis=cex.var) 
        }
        abline(h=0, lty=2)                                                 
        points(idx, coefs, pch=pch.pts, cex=cex.pts, col=col.pts)
        if (CI==2){
          segments (idx, est1[,1], idx, est1[,2], lwd=2, col=col.pts)     
          segments (idx, est2[,1], idx, est2[,2], lwd=1, col=col.pts)
        }
        else if (CI==1) {
          segments (idx, est1[,1], idx, est1[,2], lwd=1, col=col.pts)     
        }
      }
      else{
        idx <- idx + offset
        points(idx, coefs, pch=pch.pts, cex=cex.pts, col=col.pts)
        if (CI==2){
          segments (idx, est1[,1], idx, est1[,2], lwd=2, col=col.pts)     
          segments (idx, est2[,1], idx, est2[,2], lwd=1, col=col.pts)
        }
        else if (CI==1) {
          segments (idx, est1[,1], idx, est1[,2], lwd=1, col=col.pts)     
        }
      }
    }   
  }
  else{
    if (vertical){
      mar[2] <- max(min.mar[2], trunc(mar[2] + maxchar/10))  + 0.1
      par(mar=mar)
      plot(c(coefs.l, coefs.h), c(idx+k,idx-k), type="n",                                     
          axes=F, main="", xlab=xlab, ylab=ylab,...)
#      if (v.axis){
#          axis(2, n.x:1, varnames[n.x:1], las=var.las, tck=FALSE, 
#              lty=0, cex.axis=cex.var) 
#      }
    }
    else{ # horizontal
      mar[1] <- max(min.mar[1], trunc(mar[1] + maxchar/10)) + 0.1
      par(mar=mar)
      plot(c(idx+k,idx-k), c(coefs.l, coefs.h), type="n", axes=F, 
        main=main, xlab=xlab, ylab=ylab,...)                                                  
      #if (h.axis){
#          axis(1, 1:n.x, varnames[1:n.x], las=var.las, tck=FALSE, 
#              lty=0, cex.axis=cex.var) 
#      }
    }
  } 
  #on.exit(par(old.par))  
}


#' @rdname coefplot
#' @export
setMethod("coefplot", signature(object = "numeric"),
  function(object, ...)
{
  coefplot.default(object, ...)
} 
)


#' @rdname coefplot
#' @export
setMethod("coefplot", signature(object = "lm"), 
    function(object, varnames=NULL, intercept=FALSE, ...)
    {
    # collect informations
    coefs <- summary(object)$coef[,1]
    sds <- summary(object)$coef[,2]
    ifelse (is.null(varnames), varnames <- names(coefs),
            varnames <- varnames)
    if (length(varnames)!= length(names(coefs))){
      stop(message="the length of varnames does not equal the length of predictors.  
      Note: varnames must include a name for constant/intercept")
    }
    chk.int <- attr(terms(object), "intercep")
    if(chk.int & intercept | !chk.int & intercept | !chk.int & !intercept){
      intercept <- TRUE
      coefs <- coefs
      sds <- sds
      varnames <- varnames
    } else if(chk.int & !intercept){
      coefs <- coefs[-1]
      sds <- sds[-1]
      varnames <- varnames[-1]
    }    
    
    
    # plotting
    coefplot(coefs, sds, 
        varnames=varnames, ...)
    }
)
           

#' @rdname coefplot
#' @export
setMethod("coefplot", signature(object = "glm"),
    function(object, varnames=NULL, intercept=FALSE,...)
    {
    # collect informations
    coefs <- summary(object)$coef[,1]
    sds <- summary(object)$coef[,2]
    ifelse (is.null(varnames), varnames <- names(coefs),
            varnames <- varnames)
    if (length(varnames)!= length(names(coefs))){
      stop(message="the length of varnames does not equal the length of predictors.  
      Note: varnames must include a name for constant/intercept")
    }    
    chk.int <- attr(terms(object), "intercep")
    if(chk.int & intercept | !chk.int & intercept | !chk.int & !intercept){
      intercept <- TRUE
      coefs <- coefs
      sds <- sds
      varnames <- varnames
    } else if(chk.int & !intercept){
      coefs <- coefs[-1]
      sds <- sds[-1]
      varnames <- varnames[-1]
    }    
    
    # plotting
    coefplot(coefs, sds, 
        varnames=varnames, ...)
    }                                                                         
)


#' @rdname coefplot
#' @export
setMethod("coefplot", signature(object = "bugs"),
    function(object, var.idx=NULL, varnames=NULL, 
            CI=1, vertical=TRUE,
            v.axis=TRUE, h.axis=TRUE, 
            cex.var=0.8, cex.pts=0.9, 
            col.pts=1, pch.pts=20, var.las=2, 
            main=NULL, xlab=NULL, ylab=NULL, 
            plot=TRUE, add=FALSE, offset=.1,
            mar=c(1,3,5.1,2), ...)
{  
    
  if (is.null(var.idx)){
    var.idx <- 1:length(object$summary[,"50%"])
  }
  n.x <- length(var.idx)
  idx <- 1:n.x
  
  coefs <- object$summary[,"50%"][var.idx]
  if (is.null(varnames)){
    varnames <- names(coefs)     
  }
  
  if (is.null(main)){main <- "Regression Estimates"}
  if (is.null(xlab)){xlab <- ""}
  if (is.null(ylab)){ylab <- ""}
  
  min.mar <- par('mar')  
  par(mar=mar)
  
  
  maxchar <- max(sapply(varnames, nchar))

  
  k <- 1/n.x
  
  if (CI==1){
    CI50.h <- object$summary[,"75%"][var.idx]
    CI50.l <- object$summary[,"25%"][var.idx]
    CI50 <- cbind(CI50.l, CI50.h)
    if (vertical){
      mar[2] <- min(min.mar[2], trunc(mar[2] + maxchar/10)) + 0.1
      par(mar=mar)
      if(add){
        segments (CI50[,1], idx+offset, CI50[,2], idx+offset, lwd=1, col=col.pts)     
        points(coefs, idx+offset, pch=20, cex=cex.pts, col=col.pts)
      }
      else{
        plot(c(CI50[,1],CI50[,2]), c(idx+k,idx-k), type="n", 
          axes=F, main=main, xlab=xlab, ylab=ylab, ...) 
        if(plot){
          if (h.axis){
            axis(3)
          }
          if (v.axis){
          axis(2, n.x:1, varnames[n.x:1], las=var.las, tck=FALSE, 
              lty=0, cex.axis=cex.var)  
          }
          abline(v=0, lty=2)                                                 
          segments (CI50[,1], idx, CI50[,2], idx, lwd=1, col=col.pts)     
          points(coefs, idx, pch=20, cex=cex.pts, col=col.pts)  
        }
      }
    }
    else {
      mar[1] <- min(min.mar[1], trunc(mar[1] + maxchar/10))  + 0.1
      par(mar=mar)
      if(add){
          segments (idx+offset, CI50[,1], idx+offset, CI50[,2], lwd=1, col=col.pts)     
          points(idx+offset, coefs, pch=20, cex=cex.pts, col=col.pts)
      }
      else{
        plot(c(idx+k,idx-k), c(CI50[,1],CI50[,2]), type="n",                                     
            axes=F, main=main, xlab=xlab, ylab=ylab,...) 
        if(plot){
          if (v.axis){
            axis(2)
          }
          if (h.axis){
            axis(1, n.x:1, varnames[n.x:1], las=var.las, tck=FALSE, 
              lty=0, cex.axis=cex.var)  
          }
          
          abline(h=0, lty=2)                                                 
          segments (idx, CI50[,1], idx, CI50[,2], lwd=1, col=col.pts)     
          points(idx, coefs, pch=20, cex=cex.pts, col=col.pts)
        }
      }
    }   
  }
  
  if (CI==2){
    CI50.h <- object$summary[,"75%"][var.idx]
    CI50.l <- object$summary[,"25%"][var.idx]
    CI95.h <- object$summary[,"97.5%"][var.idx]
    CI95.l <- object$summary[,"2.5%"][var.idx]
    CI50 <- cbind(CI50.l, CI50.h)
    CI95 <- cbind(CI95.l, CI95.h)
    if (vertical){
      mar[2] <- min(min.mar[2], trunc(mar[2] + maxchar/10)) + 0.1
      par(mar=mar)
      if(add){
        segments (CI50[,1], idx+offset, CI50[,2], idx+offset, lwd=2, col=col.pts) 
        segments (CI95[,1], idx+offset, CI95[,2], idx+offset, lwd=1, col=col.pts)    
        points(coefs, idx+offset, pch=20, cex=cex.pts, col=col.pts)
      }
      else{
        plot(c(CI95[,1],CI95[,2]), c(idx+k,idx-k), type="n",                                     
          axes=F, main=main, xlab=xlab, ylab=ylab,...) 
        if(plot){
          if (h.axis){
            axis(3)
          }
          if (v.axis){
            axis(2, n.x:1, varnames[n.x:1], las=var.las, tck=FALSE, 
              lty=0, cex.axis=cex.var)  
          }
          abline(v=0, lty=2)                                                 
          segments (CI50[,1], idx, CI50[,2], idx, lwd=2, col=col.pts) 
          segments (CI95[,1], idx, CI95[,2], idx, lwd=1, col=col.pts)    
          points(coefs, idx, pch=20, cex=cex.pts, col=col.pts)
        }
      }
    }
    else {
      mar[1] <- min(min.mar[1], trunc(mar[1] + maxchar/10)) + 0.1
      par(mar=mar)
      if(add){
        segments (idx+offset, CI50[,1], idx+offset, CI50[,2], lwd=2, col=col.pts)
        segments (idx+offset, CI95[,1], idx+offset, CI95[,2], lwd=1, col=col.pts)         
        points(idx+offset, coefs, pch=20, cex=cex.pts, col=col.pts)        
      }
      else{
        plot(c(idx+k,idx-k), c(CI95[,1],CI95[,2]), type="n",                                     
          axes=F, main=main, xlab=xlab, ylab=ylab,...) 
        if(plot){
          if (v.axis){
            axis(2)
          }
          if (h.axis){
            axis(1, n.x:1, varnames[n.x:1], las=var.las, tck=FALSE, 
              lty=0, cex.axis=cex.var)  
          }
          abline(h=0, lty=2)                                                 
          segments (idx, CI50[,1], idx, CI50[,2], lwd=2, col=col.pts)
          segments (idx, CI95[,1], idx, CI95[,2], lwd=1, col=col.pts)         
          points(idx, coefs, pch=20, cex=cex.pts, col=col.pts)
        }
      }
    }
  }
}
)
    

#' @rdname coefplot
#' @export
setMethod("coefplot", signature(object = "polr"), 
    function(object, varnames=NULL,...)
    {
    # collect informations
    coefs <- summary(object)$coef[,1]
    sds <- summary(object)$coef[,2]
    ifelse(is.null(varnames), varnames <- names(coefs), 
        varnames <- varnames)

    # plotting
    coefplot(coefs, sds, 
        varnames=varnames, ...)
    }
)  
