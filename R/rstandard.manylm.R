#####################################################################################
# returns the standardized residuals of a multivariate linear model			#
# however standardization is NOT done multivariate, but univariate fashion 		#
#####################################################################################


rstandard.manylm <- 
function (model, infl = manylm.influence(model, do.coef = FALSE), 
    sd = sqrt(deviance(model)/df.residual(model)), ...) {

    wt.res 	<- infl$wt.res
    n 	<- NROW(wt.res)
    n.vars 	<- NCOL(wt.res)
    sD 	<- matrix(rep(sd, each=n), nrow=n, ncol=n.vars)
    hat	<- matrix(rep(infl$hat,times=n.vars), nrow=n, ncol=n.vars)
    res 	<- wt.res /(sD * sqrt(1 - hat))
    res[is.infinite(res)] <- NaN
    res

}
