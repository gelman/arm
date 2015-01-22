#' Correlation Plot
#'
#' Function for making a correlation plot starting from a data matrix.
#'
#' The function adapts the R function for Figure 8 in Tian Zheng,
#' Matthew Salganik, and Andrew Gelman, 2006, "How many people do you
#' know in prison?: using overdispersion in count data to estimate
#' social structure in networks", Journal of the American Statistical
#' Association, Vol.101, N0. 474: p.409-23.
#'
#' @param data a data matrix.
#' @param varnames variable names of the data matrix, if not provided
#'   use default variable name.
#' @param abs if TRUE, transform all correlation values into positive
#'   values, default=TRUE.
#' @param cutpts a vector of cutting points for color legend, default
#'   is NULL.  The function will decide the cutting points if cutpts is
#'   not assigned.
#' @param details show more than one digits correlaton values. Default
#'   is TRUE. FALSE is suggested to get readable output.
#' @param n.col.legend number of legend for the color thermometer.
#' @param cex.col font size of the color thermometer.
#' @param cex.var font size of the variable names.
#' @param digits number of digits shown in the text of the color theromoeter.
#' @param color color of the plot, default is FALSE, which uses gray
#'   scale.
#' @return A correlation plot.
#' @references Tian Zheng, Matthew Salganik, and Andrew Gelman, 2006,
#'   "How many people do you know in prison?: using overdispersion in
#'   count data to estimate social structure in networks", Journal of
#'   the American Statistical Association, Vol.101, N0. 474: p.409-23
#' @author
#'   Tian Zheng \email{tzheng@@stat.columbia.edu};
#'   Yu-Sung Su \email{suyusung@@tsinghua.edu.cn}
#' @seealso \code{\link[stats]{cor}}, \code{\link[graphics]{par}}
#' @export
#' @examples
#' old.par <- par(no.readonly = TRUE)
#' 
#' x1 <- rnorm(1000,50,2) 
#' x2 <- rbinom(1000,1,prob=0.63) 
#' x3 <- rpois(1000, 2) 
#' x4 <- runif(1000,40,100) 
#' x5 <- rnorm(1000,100,30)
#' x6 <- rbeta(1000,2,2) 
#' x7 <- rpois(1000,10) 
#' x8 <- rbinom(1000,1,prob=0.4) 
#' x9 <- rbeta(1000,5,4) 
#' x10 <- runif(1000,-10,-1) 
#' 
#' test.data <- data.matrix(cbind(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10))
#' test.names <- c("a short name01","a short name02","a short name03",
#'                 "a short name04","a short name05","a short name06",
#'                 "a short name07","a short name08","a short name09",
#'                 "a short name10")
#' 
#' # example 1
#' corrplot(test.data)
#' 
#' # example 2
#' corrplot(test.data,test.names, abs=FALSE, n.col.legend=7)
#' corrplot(test.data,test.names, abs=TRUE, n.col.legend=7)
#' 
#' # example 3
#' data(lalonde)
#' corrplot(lalonde, details=FALSE, color=TRUE)
#' corrplot(lalonde, cutpts=c(0,0.25,0.5,0.75), color=TRUE, digits=2)
#' 
#' par(old.par)
corrplot <- function(data, varnames=NULL, cutpts=NULL, abs=TRUE, details=TRUE, 
                     n.col.legend=5, cex.col=0.7, cex.var=0.9, digits=1, 
                     color=FALSE)
{
    
    # some check!
    if (is.matrix(data)|is.data.frame(data)){
    }
    else {
        stop ("Data must be a matrix or a data frame!")
    }
    if (sum(sapply(data, FUN=is.character))>0)
        stop ("Data contains non-numeric variables!")
    if (n.col.legend > 8)
        stop ("Suggestion: More than 8 levels of colors is difficult to read!")

    
    
    # prepare correlation matrix
    if (abs){
        z.plot <- abs(cor(data, data, use="pairwise.complete.obs"))
    }
    else{
        z.plot <- cor(data, data, use="pairwise.complete.obs")
    }
    
    if (is.null(varnames)){
        z.names <- dimnames(data)[[2]]
    }
    else{
        z.names <- varnames
    }

    triangleplot(x=z.plot, y=z.names, cutpts=cutpts, details=details, 
                n.col.legend=n.col.legend, 
                cex.col=cex.col, cex.var=cex.var, 
                digits=digits, color=color)
}
