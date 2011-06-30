################################################################################
# manyglm: Fit a multivariate generalized linear model for high dimensional data 
# the (default) methods coef, residuals, fitted values can be used             
###############################################################################

manyglm <- function (formula, family="negative.binomial", data=NULL, subset=NULL, na.action=options("na.action"), phi.method = "ML", model = FALSE, x = TRUE, y = TRUE, qr = TRUE, cor.type= "I", shrink.param=NULL, tol=1.0e-10, show.coef=FALSE, show.fitted=FALSE, show.residuals=FALSE, ... ) {

# tasmX <- as.matrix(tasmX, "numeric")  

# get family and link
#fam <- call$family
#family.char <- as.character(fam)
if ( is.character(family) ) {
    if (substr(family,1,1) == "g") {
       familynum <- 0 # gaussian
       familyname <- "gaussian"
    }   
    if (substr(family,1,1) == "p") {
       familynum <- 1   #poisson
       familyname <- "poisson"
    }   
    else if (substr(family,1,1) == "n") {
       familynum <- 2  #neg bin
       familyname <- "negative.binomial"
    }   
    else if (substr(family,1,1) == "b") {
       familynum <- 3  #logistic/binomial
       familyname <- "binomial"
    }   
    else stop (paste("'family'", family, "not recognized"))
}
else stop ("Please specify a family function with a character string. manyglm supports the following members of the exponential family: 'gaussian', 'poisson', 'binomial', 'negative.binomial' distributions.") 

# Obtain regression parameters
if (familyname=="gaussian") { 
    return(manylm(formula, data=data, subset=subset, na.action=na.action, model=model, x=x, y =y, qr=qr, cor.type=cor.type, shrink.param=shrink.param, tol=tol, ...))
}
else {
    ret.x <- x
    ret.y <- y
    ret.qr <- qr
    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)
    m  <- match(c("formula", "data", "subset", "na.action"), names(mf), 0)
    mf <- mf[c(1, m)]

    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())    # Obtain the model.frame. 
    mt <-  attr(mf, "terms")  # Obtain the model terms.

    abundances <- as.matrix(model.response(mf, "numeric"))
    if (any(is.na(abundances)) & is.null(na.action))
       stop("There are NA values in the response. An 'na.action' is necessary.")
    if(any(is.na(abundances)) & any(as.character(na.action)=="na.pass")) 
       stop("There are NA values in the response. 'na.action=na.pass' cannot be used.")

    N <- NROW(abundances)     # eg number of sites
    p <- NCOL(abundances)     # number of organism types
    labAbund<-colnames(abundances)
    if (p==0) stop("A model without response cannot be fitted.")    
    else if (p==1 & is.null(labAbund)) 
        labAbund <- deparse(attr( mt,"variables")[[2]], width.cutoff = 500) 
    else if (p>1 & is.null(labAbund)) 
        labAbund <- paste(deparse(attr( mt,"variables")[[2]],width.cutoff = 500), 1:p,sep = "") 
   labObs <- rownames(abundances)	

    Y <- abundances
    if ( familyname == "binomial" ) { # binomial distribution
       if(!is.null(labAbund) & all(substr(labAbund, 1,4) %in% c("succ", "fail")))
       {
          p <- p/2 
	  labAbund <- labAbund[substr(labAbund, 1,4) == "succ"] 	  
	  if(length(labAbund)!=p)
	     stop("for each variable of the response a column of successes and ",             "a column of failures \nmust be provided ", 
	     "(marked by 'succ' and 'fail', see 'binstruc')")
	}     
       Y[abundances>0] <- 1
    }   


    ##################### BEGIN Estimation ###################
    # Obtain the Designmatrix.
    X <- model.matrix(mt, mf)
    tX <- t(X)
    tX <- t(tX[!duplicated(tX),])  # remove duplicated col    
    if ( nrow(tX) == nrow(X) ) X <- tX
    else X <- t(tX) # for a single-column matrix 

    if (any(is.na(X)) & is.null(na.action)) 
        stop("There are NA values in the independent variables. An 'na.action' is necessary.")
    if(any(is.na(X)) & any(as.character(na.action)=="na.pass")) 
        stop("There are NA values in the independent variables. 'na.action=na.pass' cannot be used.")

    if (N!=nrow(X))
       stop("dimensions of response and independent variables in 'formula' do not match")

    q <- ncol(X)
    labX<-colnames(X)   
    # get rank
    qrx  <- qr(X)
    rank <- qrx$rank

# the following values need to be converted to integer types 
    if (phi.method == "ML") methodnum <- 0
    else if (phi.method == "Chi2") methodnum <- 1 
    else stop("Check the param estimation method name.")  

    ######### call Glm Fit Rcpp #########
    modelParam <- list(tol=tol, regression=familynum, estimation=methodnum, stablizer=FALSE)
    z <- .Call("RtoGlm", modelParam, Y, X, PACKAGE="mvabund")

# new codes added here to deal with w 
# New codes added for estimating ridge parameter 
    if (cor.type=="shrink") {      
        if (is.null(shrink.param)) {
	    tX <- matrix(1, NROW(X), 1)
            shrink.param <- ridgeParamEst(dat=z$residuals, X=tX, only.ridge=TRUE, tol=tol)$ridgeParameter
	}    
        # to simplify later computation
        if(shrink.param == 0) cor.type <- "I"
        else if(shrink.param == 1) cor.type <- "R"
        else if (abs(shrink.param)>1)
             stop("the absolute 'shrink.param' should be between 0 and 1")
    }
    else if (cor.type == "I") shrink.param <- 0 
    else if (cor.type == "R") {
         if (nrow(abundances)<=ncol(abundances))
             stop("An unstructured correlation matrix should only be used if N>>number of variables.") 
         if(nrow(abundances) < 2*ncol(abundances)) 
             warning("the calculated p-values might be unreliable as the number of cases 
                 is not much larger than the number of variables") 
         shrink.param <- 1
    }

    # parameters that are needed for tests in anova.manylm and summary.manylm.
    dimnames(z$fitted.values) <- list(labObs, labAbund)
    dimnames(z$coefficients) <- list(colnames(X), labAbund)
    dimnames(z$linear.predictor) <- list(labObs, labAbund)    
    dimnames(z$residuals) <- list(labObs, labAbund)
    dimnames(z$sqrt.1_Hii) <- list(labObs, labAbund)
    dimnames(z$var.estimator) <- list(labObs, labAbund)
    dimnames(z$sqrt.weight) <- list(labObs, labAbund)
    names(z$phi) <- labAbund
    names(z$two.loglike) <- labAbund
    names(z$deviance) <- labAbund
    names(z$aic) <- labAbund
    names(z$iter) <- labAbund

    z$tol <- tol
    z$prior.weight <- NULL
    z$AICsum <- sum(z$aic)
    z$family    <- familyname
    z$phi.method<- phi.method
    z$cor.type  <- cor.type
    z$shrink.param  <- shrink.param
    z$call      <- cl
    z$terms     <- mt
    z$formula   <- formula
    z$rank      <- rank
    z$xlevels   <- .getXlevels(mt, mf)
    z$df.residual <- N-rank
    if (model) z$model <- mf
    if (ret.x) z$x <- X
    if (ret.y) z$y <- abundances
    if (ret.qr) z$qr <- qrx 
    z$show.coef <- show.coef
    z$show.fitted <- show.fitted
    z$show.residuals <- show.residuals
 
    class(z) <- c("manyglm", "mglm")
    return(z)
  }
}

 
