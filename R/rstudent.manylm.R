#####################################################################################
# returns the studentized residuals of a multivariate linear model			#
# however standardization is NOT done multivariate, but univariate fashion 		#
#####################################################################################

rstudent.manylm <- 
function (model, infl = manylm.influence(model, do.coef = FALSE), 
    res = infl$wt.res, ...) {

    res <- res/(infl$sigma * sqrt(1 - infl$hat))
    res[is.infinite(res)] <- NaN
    res
}
