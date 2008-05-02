#R function for the logistic function
invlogit <- function (x) {
    1/(1+exp(-x))
}
