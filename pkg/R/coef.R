setMethod("coef", signature(object = "sim"),
    function(object)
    {
    ans <- object@coef
    return(ans)
    }
)

#setMethod("coef", signature(object = "sim.mer"),
#    function(object)
#    {
#    fixEffects <- object@fixef
#    ranEffects <- object@
#    
#    
#    ans <- object@coef
#    }
#)
