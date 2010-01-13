################################################################################
## anova.manylm conducts a hypotheses test comparing multivariate linear models 
################################################################################

anova.manylm <- function(object, ..., nBoot=object$nBoot,resample=object$resample, test=object$test, cor.type=object$cor.type, shrink.param=object$shrink.param, p.uni="none", studentize=TRUE, calc.rss = FALSE, tol=1.0e-10 ) 
{
    # ld.perm=TRUE load bootID from filename
    # for debug use only
    ld.perm = FALSE
    filename = NULL

    if(!any(class(object)=="manylm"))
       stop("The function 'anova.manylm' can only be used for a manylm object.")

    nRows <- nrow(object$y)
    nVars <- ncol(object$y)
    nParam <- ncol(object$x)
    Y <- matrix(as.integer(object$y), nrow=nRows, ncol=nVars) 

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
    } else if(substr(p.uni,1,1) == "s"){
       pu <- 3
       calc.pj <- adjust.pj <- TRUE
    } else
       stop("'p.uni' not defined. Choose one of 'single', 'adjusted', 'unadjusted', 'none'.")

    dimnam.a <- dimnames(object$y)[[2]]

    objects <- list(object, ...)
    nModels = length(objects)
    # ANOVA
    if (nModels>1) { 
        resdf   <- as.numeric(sapply(objects, function(x) x$df.residual))
        ####### check input arguments #######
        # check the order of models, so that each model is tested against the next smaller one
        ord <- order(resdf)
        objects <- objects[ord]
        # Test if one is nested in the other and construct the full matrix and XvarIn 
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
        #XvarIn <- matrix(as.integer(XvarIn), nrow=nModels, ncol=nParam)

        # the following values need to be converted to integer types 
        if (cor.type == "I") corr <- 1
        else if (cor.type == "R") corr <- 0
        else if (cor.type == "shrink") corr <- 2
        else stop("No such correlation type.") 

        # if (resample=="case") resam <- 0
        # To exclude case resampling
	if (resample=="case") stop("Sorry, case resampling is not yet available.")
        else if (resample == "residual") resam <- 1
        else if (resample == "score") resam <- 2
        else if (resample == "perm.resid") resam <- 3
        else stop("No such resampling method.") 
 
        if (test=="LR") testype <- 0
        else if (test == "F") testype <- 1 
        else stop("No such test method.") 
 
      # estimate ridge parameter if cor.type is not "shrink" when fitting the model
        w <- object$weights
	if (is.null(w)){
	### Fit the multivariate LM.
	   w  <- rep(1, times=nRows)
	}
	else {
	   if (!is.numeric(w))  stop("'weights' must be a numeric vector")
	   if (any(w < 0)) stop("negative 'weights' not allowed")
        }

        if (cor.type=="shrink" | is.null(shrink.param)) {
           shrink.param <- ridgeParamEst(dat=Y, X=X, weights=w, only.ridge=TRUE, doPlot=FALSE, tol=tol)$ridgeParameter
           if(shrink.param == 0) cor.type <- "I"           
           if (abs(shrink.param)>1)
              stop("the absolute 'shrink.param' should be between 0 and 1")
         }
        else if (cor.type == "I") 
           shrink.param <- 1 
        else if (cor.type == "R")
           shrink.param <- 0
        else
	    stop("no such correlation option")

        if (ld.perm) {
           if (is.null(filename)) {
               paths <- .find.package("mvabund")
              if (resample == "score")
                 filename <- file.path(paths, "data", "scores.dat")
              else if (resample == "perm.resid") 
                 filename <- file.path(paths, "data", "permID.dat")
              else
                 filename <- file.path(paths, "data", "bootID.dat")
           }
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

        ######## call resampTest Rcpp #########
        val <- .Call("RtoAnovaCpp", params, Y, X, XvarIn, bootID, PACKAGE="mvabund")

        ######## collect ANOVA results ######## 
        
        if(calc.rss) {
           objects <- list()
           for (i in 1:nModels) {
               Xi <- X[, which(XvarIn[i, ]>0)]
               objects[[i]] <- lm(Y~Xi)
           }
           RSS <- matrix(unlist(lapply(objects, deviance)),nrow=nModels,ncol=nVars,byrow=TRUE)
           if (is.null(dimnam.a)) dimnam.a <- paste("abund", 1:nVars)
           dimnames(RSS) <- list(paste("Model", 1:nModels) , dimnam.a)
           Diff <- matrix(unlist(lapply(1:nVars, function(x) diff(RSS[ord,x]) )),           nrow=nModels-1, ncol=nVars,byrow=TRUE )
           dimnames(Diff) <- list(paste("Model", 1:(nModels-1), 2:nModels), dimnam.a)
           attr(RSS, "title") <- "\nResidual Sum of Squares\n"
           attr(Diff,"title") <- "\nDiff. Sum of Squares\n"
        }
        else  RSS <- Diff <- NULL
 
         anova <- list()
         # Outputs passed from inputs
         anova$p.uni <- p.uni
         anova$test  <- test
         anova$cor.type <- cor.type
         anova$resample <- resample
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
         Xnames <- lapply(objects, function(x) paste(deparse(formula(x), width.cutoff=500), collapse = "\n")) 
         topnote <- paste("Model ", format(1:nModels), ": ", Xnames, sep = "", collapse = "\n")
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
         dimnames(anova$table) <- list(paste("Model", 1:nModels), c("Res.Df", "Df", testname, pname))
         # make several univariate tables 
         attr(anova$uni.test, "title") <- attr(anova$uni.p, "title") <- "\nUnivariate Tests\nTest statistics:\n"
         if (is.null(dimnam.a)) dimnam.a <- paste("abund", 1:nVars)
         dimnames(anova$uni.p) <- dimnames(anova$uni.test) <- list(paste("Model",(1:nModels)), dimnam.a)

         class(anova) <- "anova.manylm"
         return(anova)
     }
  else {
  #  return the summary multstat result
     smry <- summary.manylm(object, nBoot=nBoot,resample=resample, test=test, cor.type=cor.type, shrink.param=shrink.param, p.uni=p.uni, studentize=studentize, tol=tol) 

         anova <- list()
         # Outputs passed from inputs
         anova$p.uni <- p.uni
         anova$test  <- test
         anova$cor.type <- cor.type
         anova$resample <- resample
         anova$shrink.param <- shrink.param
         anova$nBoot <- nBoot 
         # parameter
         anova$calc.rss <- FALSE 
         anova$n.bootsdone <- smry$nSamp
         anova$n.iter.sign <- nBoot - smry$nSamp
         anova$one <- TRUE
         # model fit
         anova$RSS   <- NULL 
         # test statistics
         rank  <- object$rank
         df.residual <- object$df.residual
         if (rank==0) { # i.e. y~0 --> a more general model is not possible
             anova$table<- as.data.frame(matrix(NA, 1, 1))
             anova$table[1,1] <- df.residual
             dimnames(anova$table) <- list( "Residuals" , c("Df") )
             attr(anova$table, "heading") <- c(title, topnote)
             if (calc.pj){
                anova$uni.test <- anova$uni.p <- matrix(ncol=nVars,nrow=0)
                colnames(anova$uni.test) <- colnames(anova$uni.p) <-                         colnames(object$y)
             }
             anova$test <- NULL 
             class(anova) <- "anova.manylm"
             return(anova)
        }
         anova$table <- list() 
         anova$table <- data.frame(c(smry$statistic[3], object$df), c(smry$statistic[1], NA), c(smry$statistic[2], NA))
         # display formats        
        topnote <- paste("Response: ", deparse(formula(object)[[2]], width.cutoff = 500), collapse = "\n")
        title <- "Analysis of Variance Table \n" 
        attr(anova$table, "heading") <- c(title, topnote)
        obj.terms <- object$terms
        nmeffects <- c(attr(obj.terms, "term.labels"))
        p1 <- 1:rank
        asgn <- object$assign[object$qr$pivot][p1]
        indVarNames <- nmeffects[ unique(asgn)]

        attr(anova$table, "title") <-
        "\nOverall test for all response variables\nTest statistics:\n" 
        if (!is.null(test)) {
           testname <- paste(test,"value")
           pname    <- paste("Pr(>",test,")", sep="")
        } else {
           testname <- "no test"
           pname    <- ""
        }
        dimnames(anova$table)[[1]] <- c(indVarNames, "Residuals")
        dimnames(anova$table)[[2]][1] <- c("Df")
        dimnames(anova$table)[[2]][2:3] <- c(testname, pname)
        if (calc.pj){ 
            sp <- matrix(smry$statistic.j[, 2], ncol=nVars, nrow=1)
            st <- matrix(smry$statistic.j[, 1], ncol=nVars, nrow=1)
            anova$uni.p <- sp
            anova$uni.test <- st
            attr(anova$uni.p, "title") <- "\nUnivariate Tests\nP values:\n"
            attr(anova$uni.test, "title") <- "\nUnivariate Tests\nTest statistics:\n"
            dimnames(anova$uni.p) <- dimnames(anova$uni.test) <- list(indVarNames, dimnam.a)    
        }
        else {
            anova$uni.p <- NULL
            anova$uni.test <- NULL
        }
       
        class(anova) <- "anova.manylm"
        return(anova)
  }
}
