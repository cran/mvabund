################################################################################
# PLOT.MANYLM, PLOT.MANYGLM:                                                   #
# Plot for evaluation of goodness of fit for lm.mvabund objects                #
################################################################################

plot.manyglm <- plot.manylm  <- function(x, 
					which = 1:4, 
					caption = c("Residuals vs Fitted","Normal Q-Q", "Scale-Location", "Cook\'s distance"), 
					overlay=TRUE,
    					n.vars=12, 
					var.subset=NULL, 
					sub.caption = NULL, 
					studentized=TRUE,
					... ) 


{	 
	#Pass argument to appropriate function.
	isGlm <- inherits(x,"manyglm")

	if (isGlm) {

		default.plot.manyglm(x, 
				which = which,
				caption = caption, 
				overlay=overlay,
    				n.vars=n.vars, 
				var.subset=var.subset, 
				sub.caption = sub.caption, 
				studentized=studentized,
				... )
		
		return(invisible())
	} else {

		default.plot.manylm(x, 
				which = which, 
				caption = caption, 
				overlay=overlay,
    				n.vars=n.vars, 
				var.subset=var.subset, 
				sub.caption = sub.caption, 
				studentized=studentized,
				... )

		return(invisible())

	}

}
 
