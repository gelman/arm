#' @rdname sim
#' @export
setMethod("coef", signature(object = "sim"),
    function(object)
    {
    ans <- object@coef
    return(ans)
    }
)

#' @rdname sim
#' @export
setMethod("coef", signature(object = "sim.polr"),
    function(object, slot=c("ALL", "coef", "zeta"))
    {
    slot <- match.arg(slot)
    if(slot=="coef"){
      ans <- object@coef
    } else if(slot=="zeta"){
      ans <- object@zeta
    } else {
      ans <- cbind(object@zeta, object@coef)
    }
    return(ans)
    }
)

#' @rdname sim
#' @export
setMethod("coef", signature(object = "sim.merMod"),
    function(object)
    {
    fef <- object@fixef
    ref <- object@ranef
    ans <- list("fixef" = fef, "ranef" = ref)
    return(ans)
    }
)


#' @rdname sim
#' @export
fixef.sim.merMod <- function(object, ...) {
  object@fixef
}

#' @rdname sim
#' @export
ranef.sim.merMod <- function(object, ...) {
  object@ranef
}
