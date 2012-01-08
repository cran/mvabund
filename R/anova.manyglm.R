###############################################################################
# R user interface to anova test for comparing multivariate linear models 
# Author: Yi Wang (yi dot wang at computer dot org)
# 11-Nov-2011
###############################################################################

anova.manyglm <- function(object, ..., resamp="montecarlo", test="LR", p.uni="none", nBoot=1000, cor.type=object$cor.type, ld.perm=FALSE, filename=NULL ) 
{
    if (cor.type!="I" & test=="LR") {
        warning("The likelihood ratio test can only be used if correlation matrix of the abundances is is assumed to be the Identity matrix. The Wald Test will be used.")
        test <- "wald"
    }

    if (any(class(object) == "manylm")) {
        if ( test == "LR" ) 
	    return(anova.manylm(object, ..., resamp=resamp, test="LR", p.uni=p.uni, nBoot=nBoot, cor.type=cor.type, shrink.param=object$shrink.param, tol=tol, ld.perm=ld.perm, filename=filename))
        else {
	    warning("For an manylm object, only the likelihood ratio test and F test are supported. So the test option is changed to `'F''. ")
	    return(anova.manylm(object, resamp=resamp, test="F", p.uni=p.uni, nBoot=nBoot, cor.type=cor.type, tol=tol, ld.perm=ld.perm, filename=filename, ... ))
	}
    }   
    else if (!any(class(object)=="manyglm"))
        stop("The function 'anova.manyglm' can only be used for a manyglm or manylm object.")

    #check if any non manylm object in ...
    objects <- list(object, ...)    
    dots <- list(...)
    ndots <- length(dots)
    if (ndots>0) {
       which <- rep(TRUE, ndots+1)
       for (i in 1:ndots) {
           if (!any(class(dots[[i]])=="manyglm")){
              objectname <- names(dots[i])
              warning(paste(objectname, "is not a manyglm object nor a valid argument- removed from input"))
              which[i+1] <- FALSE
           }
       }
       objects <- objects[which]
    }

    tol = object$tol
    nModels = length(objects)
    nRows <- nrow(object$y)
    nVars <- ncol(object$y)
    nParam <- ncol(object$x)
    dimnam.a <- dimnames(object$y)[[2]]
    if (is.null(dimnam.a)) dimnam.a <- paste("abund", 1:nVars)

    Y <- matrix(as.integer(object$y), nrow=nRows, ncol=nVars) 
    if (is.null(Y)) {
    #      mu.eta <- object$family$mu.eta
        eta <- object$linear.predictor
        Y <- object$fitted.values + object$residuals * log(eta)	     
     }

    w <- object$weights
    if (is.null(w)) w  <- rep(1, times=nRows)
    else {
        if (!is.numeric(w))  stop("'weights' must be a numeric vector")
        if (any(w < 0)) stop("negative 'weights' not allowed")
    }

    # the following values need to be converted to integer types  
    if (substr(object$family,1,1) == "p") familynum <- 1 
    else if (substr(object$family,1,1) == "n") familynum <- 2
    else if (substr(object$family,1,1)=="b") familynum <- 3

    if (object$phi.method == "ML") methodnum <- 0
    else if (object$phi.method == "Chi2") methodnum <- 1 

    if (substr(resamp,1,1)=="c") resampnum <- 0  #case
    # To exclude case resampling
    #if (resamp=="case") stop("Sorry, case resampling is not yet available.")
    else if (substr(resamp,1,1)=="r") resampnum <- 1  # residual
    else if (substr(resamp,1,1)=="s") resampnum <- 2  # score
    else if (substr(resamp,1,1) =="p") resampnum <- 3 # permuation
#    else if (substr(resamp,1,1) =="f") resampnum <- 4 # free permuation
    else if (substr(resamp,1,1) ==  "m") resampnum <- 5 # montecarlo 
    else stop("'resamp' not defined. Choose one of 'case', 'resid', 'score', 'perm.resid', 'montecarlo'")    
        
    if (substr(test,1,1) == "w") testnum <- 2 # wald
    else if (substr(test,1,1) == "s") testnum <- 3 #score
    else if (substr(test,1,1) == "L") testnum <- 4 #LR
    else stop("'test'not defined. Choose one of 'wald', 'score', 'LR' for an manyglm object.")  

    if (cor.type == "R") corrnum <- 0
    else if (cor.type == "I") corrnum <- 1
    else if (cor.type == "shrink") corrnum <- 2
    else stop("'cor.type' not defined. Choose one of 'I', 'R', 'shrink'")  

    if (ld.perm && !is.null(filename)) {
        bootID <- as.matrix(read.table(filename), nrow=nBoot, ncol=nRows)
        rep <- 1
    }
    else {
        bootID <- c(FALSE)
        rep <- 0
    }

    if(substr(p.uni,1,1) == "n"){
       pu <- 0
       calc.pj <- adjust.pj <- FALSE
    } else if(substr(p.uni,1,1) == "u"){
       pu <- 1
       calc.pj <- TRUE
       adjust.pj <- FALSE
    } else if(substr(p.uni,1,1)=="a"){
       pu <- 2
       calc.pj <- adjust.pj <- TRUE
    } else
       stop("'p.uni' not defined. Choose one of 'adjusted', 'unadjusted', 'none'.")

    # construct for param list     
    modelParam <- list(tol=tol, regression=familynum, 
                       estimation=methodnum, stablizer=0)
    # note that nboot excludes the original data set
    testParams <- list(tol=tol, nboot=nBoot-1, cor_type=corrnum, 
              test_type=testnum, resamp=resampnum, reprand=rep, punit=pu)

    # ANOVA
    if (nModels==1) {
       # test the significance of each model terms
       X <- object$x
       varseq <- object$assign
       nterms <- max(0, varseq)+1
       resdev <- resdf <- NULL
       XvarIn <- matrix(ncol=nParam, nrow=nterms, 1)
       for ( i in 0L:(nterms-2)) { # exclude object itself
           XvarIn[nterms-i, varseq>i] <- 0 # in reversed order
           ncoef <- nParam-length(varseq[varseq>i])
           resdf <- c(resdf, nRows-ncoef)
       }
       resdf <- c(resdf, object$df.residual)
#browser()       
       # get the shrinkage estimates
       tX <- matrix(1, nrow=nRows, ncol=1)
       if (corrnum==2 | resampnum==5){ # shrinkage or montecarlo bootstrap
          shrink.param <- c(rep(NA, nterms))
          if (object$cor.type == "shrink") 
       	      shrink.param[1] <- object$shrink.param	 
	  else    
              shrink.param[1] <- ridgeParamEst(dat=object$residuals, X=tX, 
	                     only.ridge=TRUE)$ridgeParam
          for ( i in 0:(nterms-2)){ # exclude object itself
              fit <- .Call("RtoGlm", modelParam, Y, X[,varseq<=i,drop=FALSE], 
	              PACKAGE="mvabund")
              shrink.param[nterms-i] <- ridgeParamEst(dat=fit$residuals, 
                      X=tX, only.ridge=TRUE)$ridgeParam # in reversed order
           }
       }   
       else if (corrnum == 0) shrink.param <- c(rep(1, nterms))
       else if (corrnum == 1) shrink.param <- c(rep(0, nterms))
#       resdev <- c(resdev, object$deviance) 
       nModels <- nterms

       ord <- (nterms-1):1
       topnote <- paste("Model:", deparse(object$call))
       tl <- attr(object$terms, "term.labels")
       tl <- c("(Intercept)", tl)
    }   
    else {
        targs <- match.call(call = sys.call(which = 1), expand.dots = FALSE)
        if (targs[[1]] == "example")
            modelnamelist <- paste("Model ", format(1:nModels))
        else    
            modelnamelist <- as.character(c(targs[[2]], targs[[3]]))

        resdf   <- as.numeric(sapply(objects, function(x) x$df.residual))
        ####### check input arguments #######
        # check the order of models, so that each model is tested against the next smaller one 
        ord <- order(resdf, decreasing=TRUE) 
        objects <- objects[ord]
        resdf <- resdf[ord]
        modelnamelist <- modelnamelist[ord]

        # get the shrinkage estimates
        if (corrnum == 2 | resampnum == 5) { # shrinkage or parametric bootstrap
	    shrink.param <- c(rep(NA,nModels))
    	    tX <- matrix(1, nrow=nRows, ncol=1)
	    for ( i in 1:nModels ) {
	        if (objects[[i]]$cor.type == "shrink") 
                    shrink.param[i] <- objects[[i]]$shrink.param
	        else shrink.param[i] <- ridgeParamEst(dat=objects[[i]]$residuals, X=tX, only.ridge=TRUE)$ridgeParam 
	    }
	}
        else if (corrnum == 0) shrink.param <- c(rep(1,nModels))
        else if (corrnum == 1) shrink.param <- c(rep(0,nModels))

        # Test if models are nested, construct the full matrix and XvarIn 
        XNull <- as.matrix(objects[[1]]$x, "numeric")
        ind <- matrix(ncol=1, nrow=nModels)
        for ( i in 2:nModels ) {
            XAlt  <- as.matrix(objects[[i]]$x, "numeric")
            Xarg  <- cbind(XAlt, XNull)
            tmp <- qr(Xarg)
            Xplus <- qr(XAlt)
            if ( tmp$rank == Xplus$rank ) {
               Beta <- qr.coef(Xplus, XNull)  # equivalent to (XAlt\XNull) in matlab 
               # The following gets the left null space of beta, ie.LT=null(t(beta));
               # note that LT is an orthogonal complement of Beta, and [Beta, LT] together forms the orthogonal basis that span the column space of XAlt
               # For some reason, it must be null(beta) instead of null(t(beta)) in R to get the same answer in matlab.
               tmp <- qr(Beta)
               set <- if(tmp$rank == 0) 1:ncol(Beta) else  - (1:tmp$rank)
               LT <- qr.Q(tmp, complete = TRUE)[, set, drop = FALSE]
               # to get the dimension of Xnull
               ind[nModels+2-i, 1] <- dim(XNull)[2]
               XNull <- cbind(XNull, XAlt%*%LT)
            } 
            else
              stop(paste(modelnamelist[i-1], "is note nested in Model", modelnamelist[i]))
        }
        # the full matrix template X, note that Xnull and Xalt are reconstructed from X and XvarIn in the resampling process
        X <- XNull
        nParam <- ind[1, 1] <- dim(X)[2] 
        XvarIn <- matrix(ncol=nParam, nrow=nModels, as.integer(0))  
        Xnames <- list()   # formula of each model
        for ( i in 1:nModels ) XvarIn[i, 1:ind[i, 1]] <- as.integer(1) 

        Xnames <- lapply(objects, function(x) paste(deparse(formula(x), 
                         width.cutoff=500), collapse = "\n")) 
        topnote <- paste(modelnamelist, ": ", Xnames, sep = "", collapse = "\n")
        tl <- modelnamelist	
        ord <- (nModels-1):1
     }
#browser()
    ######## call resampTest Rcpp #########
    val <- .Call("RtoGlmAnova", modelParam, testParams, Y, X, 
                 XvarIn, bootID, shrink.param, PACKAGE="mvabund")

    # prepare output summary
    table <- data.frame(resdf, c(NA, val$dfDiff[ord]), 
                 c(NA, val$multstat[ord]), c(NA, val$Pmultstat[ord])) 
    uni.p <- matrix(ncol=nVars,nrow=nModels) 
    uni.test <- matrix(ncol=nVars, nrow=nModels)
    uni.p[2:nModels, ] <- val$Pstatj[ord,]   
    uni.test[2:nModels, ] <- val$statj[ord,]

    anova <- list()
    # Supplied arguments
    anova$family <- object$family
    anova$p.uni <- p.uni
    anova$test  <- if (test=="LR") "Dev" else test
    anova$cor.type <- cor.type
    anova$resamp <- if (resamp=="montecarlo") "parametric" else resamp
    anova$nBoot <- nBoot 
    # estimated parameters
    anova$shrink.param <- shrink.param
    anova$n.bootsdone <- val$nSamp
    # test statistics
    anova$table <- table 
    anova$uni.p <- uni.p
    anova$uni.test <- uni.test

    ########### formal displays #########
    # Title and model formulas
    title <- if (test=="LR") "Analysis of Deviance Table\n"
             else "Analysis of Variance Table\n" 
    attr(anova$table, "heading") <- c(title, topnote) 
    attr(anova$table, "title") <- "\nMultivariate test:\n" 
    # make multivariate table 
    if (!is.null(test)) {
       testname <- anova$test
       pname    <- paste("Pr(>",anova$test,")", sep="")
    } else {
       testname <- "no test"
       pname    <- ""
    }
    
    dimnames(anova$table) <- list(tl, c("Res.Df", "Df.diff", testname, pname))
   
    # make several univariate tables 
    attr(anova$uni.test, "title") <- attr(anova$uni.p, "title") <- "\nUnivariate Tests:\n"
    dimnames(anova$uni.p) <- dimnames(anova$uni.test) <- list(tl, dimnam.a)

    class(anova) <- "anova.manyglm"
    return(anova)
}
