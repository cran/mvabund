###############################################################################
# R user interface to anova test for comparing multivariate linear models 
# Author: Yi Wang (yi dot wang at computer dot org)
# 05-Jan-2010
###############################################################################

anova.manylm <- function(object, ..., resamp="perm.resid", test="F", p.uni="none", nBoot=1000, cor.type=object$cor.type, shrink.param=object$shrink.param, studentize=TRUE, calc.rss = FALSE, tol=1.0e-10, ld.perm=FALSE, filename=NULL ) 
{
    if(!any(class(object)=="manylm"))
       stop("The function 'anova.manylm' can only be used for a manylm object.")

    #check if any non manylm object in ...
    objects <- list(object, ...)
    dots <- list(...)
    ndots <- length(dots)
    if (ndots > 0) {
       which <- rep(TRUE, ndots+1)
       for (i in 1:ndots) {
           if (!any(class(dots[[i]])=="manylm")){
              objectname <- names(dots[i])
              warning(paste(objectname, "is not a manylm object nor a valid argument- removed from input")) 
              which[i+1] <- FALSE
           }           
       }        
       objects <- objects[which]
    }

    dimnam.a <- dimnames(object$y)[[2]]
    if (is.null(dimnam.a)) dimnam.a <- paste("abund", 1:nVars)
    nRows <- nrow(object$y)
    nVars <- ncol(object$y)
    nParam <- ncol(object$x)
    nModels = length(objects)
    Y <- matrix(as.integer(object$y), nrow=nRows, ncol=nVars) 
    X <- matrix(object$x, nrow=nRows, ncol=nParam)

    # the following values need to be converted to integer types 
    w <- object$weights
    if (is.null(w)) w <- rep(1, times=nRows)
    else {
        if (!is.numeric(w))  stop("'weights' must be a numeric vector")
        if (any(w < 0)) stop("negative 'weights' not allowed")
    }

    if (substr(p.uni,1,1) == "n"){
       pu <- 0
       calc.pj <- adjust.pj <- FALSE
    } else if(substr(p.uni,1,1) == "u"){
       pu <- 1
       calc.pj <- TRUE
       adjust.pj <- FALSE
    } else if(substr(p.uni,1,1)=="a"){
       pu <- 2
       calc.pj <- adjust.pj <- TRUE
    } else if(substr(p.uni,1,1) == "s"){
       pu <- 3
       calc.pj <- adjust.pj <- TRUE
    } else
       stop("'p.uni' not defined. Choose one of 'single', 'adjusted', 'unadjusted', 'none'.")

    if (resamp=="case") resam <- 0
    # To exclude case resampling
    # if (resamp=="case") 
    #   stop("Sorry, case resampling is not yet available.")
    else if (resamp == "residual") resam <- 1
    else if (resamp == "score") resam <- 2
    else if (resamp == "perm.resid") resam <- 3
    else stop("No such resampling method.") 
 
    if (test=="LR") testype <- 0
    else if (test == "F") testype <- 1 
    else stop("No such test method.") 

    if (cor.type == "I") {
       corr <- 1
       shrink.param <- 0
    }    
    else if (cor.type == "R") { 
       corr <- 0 
       shrink.param <- 1
    }
    else if (cor.type == "shrink") corr <- 2
    else stop("No such correlation type.") 

   if (is.null(shrink.param)) {
       if ( object$cor.type=="shrink" ) 
          shrink.param=object$shrink.param
       else shrink.param <- ridgeParamEst(dat=Y, X=X, weights=w, only.ridge=TRUE, tol=tol)$ridgeParameter
       if (abs(shrink.param)>1)
          stop("the absolute 'shrink.param' should be between 0 and 1")
   }

    if (ld.perm && !is.null(filename)) {
        bootID <- as.matrix(read.table(filename), nrow=nBoot, ncol=nRows)
        rep <- 1
    }
    else {
        bootID <- c(FALSE)
        rep <- 0
    }

    if (studentize) st <- 1
    else st <- 0

    # construct for param list      
    params <- list(tol=tol, nboot=nBoot, cor_type=corr, shrink_param=shrink.param, test_type=testype, resamp=resam, reprand=rep, studentize=st, punit=pu, rsquare=0)

    # ANOVA
    if (nModels==1) {
        varseq <- object$assign
        nterms <- max(0, varseq)+1
        resdf  <- resdev <- NULL

        XvarIn <- matrix(ncol=nParam, nrow=nterms, 1)
        for ( i in 0L:(nterms-2)){ # exclude object itself
            XvarIn[nterms-i, varseq>i] <- 0 # in reversed order
            Xi <- X[, varseq<=i, drop=FALSE]
            fit <- manylm(Y~Xi)
            deviance <- as.numeric(deviance.manylm(fit))
            resdev <- c(resdev, deviance)
            resdf <- c(resdf, nRows-dim(fit$coefficients)[1])
        }
        resdf <- c(resdf, object$df.residual)
        deviance <- as.numeric(deviance.manylm(object))
        resdev <- c(resdev, deviance)

        nModels <- nterms

        ord <- nterms:1
        topnote <- paste("Model:", deparse(object$call) )
        tl <- attr(object$terms, "term.labels")
        tl <- c("Null", tl)
    }
    else {
        resdf <- as.numeric(sapply(objects, function(x) x$df.residual))
        resdev <- as.numeric(sapply(objects, function(x) deviance.manylm(x)))
        ####### check input arguments #######
        # each model is tested against the next smaller one
        ord <- order(resdf)
        objects <- objects[ord]
        # construct a list of nested models
        XNull <- as.matrix(objects[[nModels]]$x, "numeric")
        ind <- matrix(ncol=1, nrow=nModels)
        for ( i in nModels-c(2:nModels-1) ) {
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
               ind[i+1, 1] <- dim(XNull)[2]
               XNull <- cbind(XNull, XAlt%*%LT)
            } 
            else
               stop(paste("Model", ord[i+1], "is note nested in Model", ord[i]))
        }
        # the full matrix template X, note that Xnull and Xalt are reconstructed from X and XvarIn in the resampling process
        X <- XNull
        nParam <- ind[1, 1] <- dim(X)[2] 
        XvarIn <- matrix(ncol=nParam, nrow=nModels, as.integer(0))  
        Xnames <- list()   # formula of each model
        for ( i in 1:nModels ) XvarIn[i, 1:ind[i, 1]] <- as.integer(1) 

        Xnames <- lapply(objects, function(x) paste(deparse(formula(x), 
                  width.cutoff=500), collapse = "\n")) 
        topnote <- paste("Model ", format(1:nModels), ": ", 
                  Xnames, sep = "", collapse = "\n")
     
        tl <- paste("Model", 1:nModels)
    }

    ######## call resampTest Rcpp #########
    val <- .Call("RtoAnovaCpp", params, Y, X, XvarIn, bootID, PACKAGE="mvabund")

    if (calc.rss) {
        RSS <- matrix(unlist(resdev),nrow=nModels,ncol=nVars,byrow=TRUE)
        dimnames(RSS) <- list(paste("Model", 1:nModels) , dimnam.a)
        Diff <- matrix(unlist(lapply(1:nVars, function(x) diff(RSS[ord,x]) )),          nrow=nModels-1, ncol=nVars,byrow=TRUE )
        dimnames(Diff) <- list(paste("Model", 1:(nModels-1), 2:nModels), dimnam.a)
        attr(RSS, "title") <- "\nResidual Sum of Squares\n"
        attr(Diff,"title") <- "\nDiff. Sum of Squares\n"
    }
    else  RSS <- Diff <- NULL
 
    ######## collect ANOVA results ######## 
    anova <- list()
    # Outputs passed from inputs
    anova$p.uni <- p.uni
    anova$test  <- test
    anova$cor.type <- cor.type
    anova$resamp <- resamp
    anova$shrink.param <- shrink.param
    anova$nBoot <- nBoot 
    # parameter
    anova$calc.rss <- calc.rss  
    anova$n.bootsdone <- val$nSamp
    anova$n.iter.sign <- nBoot - val$nSamp
    anova$one <- FALSE
    # model fit
    anova$RSS   <- list()
    anova$RSS$RSS <- RSS
    anova$RSS$Diff <- Diff
    # test statistics
    anova$table <- data.frame(resdf[ord], c(NA, val$dfDiff), c(NA, val$multstat), c(NA, val$Pmultstat)) 
    anova$uni.p <- matrix(ncol=nVars,nrow=nModels) 
    anova$uni.test <- matrix(ncol=nVars, nrow=nModels)
    anova$uni.p[2:nModels, ] <- val$Pstatj
    anova$uni.test[2:nModels, ] <- val$statj

    ########### formal displays #########
    # Title and model formulas
    title <- "Analysis of Variance Table\n" 
    attr(anova$table, "heading") <- c(title, topnote) 
    attr(anova$table, "title") <- "\nOverall test for all response variables\nTest statistics:\n" 

    # make multivariate table 
    if (!is.null(test)) {
       testname <- paste("val(",test,")", sep="")
       pname    <- paste("Pr(>",test,")", sep="")
    } else {
       testname    <- "no test"
       pname       <- ""
    }
    dimnames(anova$table) <- list(tl, c("Res.Df", "Df.diff", testname, pname))
    # make several univariate tables 
    attr(anova$uni.test, "title") <- attr(anova$uni.p, "title") <- "\nUnivariate Tests\nTest statistics:\n"
    dimnames(anova$uni.p) <- dimnames(anova$uni.test) <- list(tl, dimnam.a)

    class(anova) <- "anova.manylm"
    return(anova)
}
