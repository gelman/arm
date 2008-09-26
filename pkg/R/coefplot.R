coefplot.default <- function(coefs, sds, 
            varnames=NULL, CI=2, 
            vertical=TRUE,
            v.axis=TRUE, h.axis=TRUE,
            cex.var=0.8, cex.pts=0.9, col.pts=1, pch.pts=20,
            var.las=2, main=NULL, xlab=NULL, ylab=NULL, mar=c(1,3,5.1,2),
            plot=TRUE, add=FALSE, epsilon=0.1,...)
{
  
     # collect informations
    if (is.list(coefs)){
      coefs <- unlist(coefs)
    }
    n.x <- length(coefs)
    idx <- seq(1, n.x)   
    
    coefs.h <- coefs + CI*sds 
    coefs.l <- coefs - CI*sds                                                          
    
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
            segments (coefs+sds, idx, coefs-sds, idx, lwd=2, col=col.pts)     
            segments (coefs+2*sds, idx, coefs-2*sds, idx, lwd=1, col=col.pts)
          }
          else{
            segments (coefs+sds, idx, coefs-sds, idx, lwd=1, col=col.pts)    
          }
        }
        else{
          idx <- idx + epsilon
          points(coefs, idx, pch=pch.pts, cex=cex.pts, col=col.pts)
          if (CI==2){
            segments (coefs+sds, idx, coefs-sds, idx, lwd=2, col=col.pts)     
            segments (coefs+2*sds, idx, coefs-2*sds, idx, lwd=1, col=col.pts)
          }
          else{
            segments (coefs+sds, idx, coefs-sds, idx, lwd=1, col=col.pts)    
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
          segments (idx, coefs+sds, idx, coefs-sds, lwd=2, col=col.pts)     
          segments (idx, coefs+2*sds, idx, coefs-2*sds, lwd=1, col=col.pts)
        }
        else if (CI==1) {
          segments (idx, coefs+sds, idx, coefs-sds, lwd=1, col=col.pts)     
        }
      }
      else{
        idx <- idx + epsilon
        points(idx, coefs, pch=pch.pts, cex=cex.pts, col=col.pts)
        if (CI==2){
          segments (idx, coefs+sds, idx, coefs-sds, lwd=2, col=col.pts)     
          segments (idx, coefs+2*sds, idx, coefs-2*sds, lwd=1, col=col.pts)
        }
        else if (CI==1) {
            segments (idx, coefs+sds, idx, coefs-sds, lwd=1, col=col.pts)     
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
}

setMethod("coefplot", signature(object = "numeric"),
  function(object, ...)
{
  coefplot.default(object, ...)
} 
)



setMethod("coefplot", signature(object = "lm"), 
    function(object, varnames=NULL, intercept=FALSE, add=FALSE, ...)
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
    if (intercept){
        coefs <- coefs
        sds <- sds
        varnames <- varnames
    }
    else{
        coefs <- coefs[-1]
        sds <- sds[-1]
        varnames <- varnames[-1]
    }
    
    
    # plotting
    coefplot(coefs, sds, 
        varnames=varnames, ...)
    }
)
           
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
    if (intercept){
        coefs <- coefs
        sds <- sds
        varnames <- varnames
    }
    else{
        coefs <- coefs[-1]
        sds <- sds[-1]
        varnames <- varnames[-1]
    }
    
    
    # plotting
    coefplot(coefs, sds, 
        varnames=varnames, ...)
    }                                                                         
)


setMethod("coefplot", signature(object = "bugs"),
    function(object, var.idx=NULL, varnames=NULL, 
            CI=1, vertical=TRUE,
            v.axis=TRUE, h.axis=TRUE, 
            cex.var=0.8, cex.pts=0.9, 
            col.pts=1, pch.pts=20, var.las=2, 
            main=NULL, xlab=NULL, ylab=NULL, 
            plot=TRUE, add=FALSE, epsilon=.1,
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
      mar[2] <- max(min.mar[2], trunc(mar[2] + maxchar/10)) + 0.1
      par(mar=mar)
      if(add){
        segments (CI50[,1], idx+epsilon, CI50[,2], idx+epsilon, lwd=1, col=col.pts)     
        points(coefs, idx+epsilon, pch=20, cex=cex.pts, col=col.pts)
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
      mar[1] <- max(min.mar[1], trunc(mar[1] + maxchar/10))  + 0.1
      par(mar=mar)
      if(add){
          segments (idx+epsilon, CI50[,1], idx+epsilon, CI50[,2], lwd=1, col=col.pts)     
          points(idx+epsilon, coefs, pch=20, cex=cex.pts, col=col.pts)
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
      mar[2] <- max(min.mar[2], trunc(mar[2] + maxchar/10)) + 0.1
      par(mar=mar)
      if(add){
        segments (CI50[,1], idx+epsilon, CI50[,2], idx+epsilon, lwd=2, col=col.pts) 
        segments (CI95[,1], idx+epsilon, CI95[,2], idx+epsilon, lwd=1, col=col.pts)    
        points(coefs, idx+epsilon, pch=20, cex=cex.pts, col=col.pts)
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
      mar[1] <- max(min.mar[1], trunc(mar[1] + maxchar/10)) + 0.1
      par(mar=mar)
      if(add){
        segments (idx+epsilon, CI50[,1], idx+epsilon, CI50[,2], lwd=2, col=col.pts)
        segments (idx+epsilon, CI95[,1], idx+epsilon, CI95[,2], lwd=1, col=col.pts)         
        points(idx+epsilon, coefs, pch=20, cex=cex.pts, col=col.pts)        
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
