################################################################################
# manyglm: Fit a multivariate generalized linear model for high dimensional data
# the (default) methods coef, residuals, fitted values can be used
###############################################################################

manyglm <- function (formula,
    family="negative.binomial",
    composition=FALSE,
    data=NULL,
    subset=NULL,
    na.action=options("na.action"),
    K=1,
    theta.method = "PHI",
    model = FALSE,
    x = TRUE,
    y = TRUE,
    qr = TRUE,
    cor.type= "I",
    shrink.param=NULL,
    tol=sqrt(.Machine$double.eps),
    maxiter=25,
    maxiter2=10,
    show.coef=FALSE, show.fitted=FALSE, show.residuals=FALSE, show.warning=FALSE,
    offset, ... ) {


# start by converting any family objects that can be handled over to character strings
if(class(family)=="family"){
  fam = family
  if(fam$family=="binomial" & fam$link=="cloglog")
    family="cloglog"
  if(fam$family=="binomial" & fam$link=="logit")
    family="binomial"
  if(fam$family=="poisson" & fam$link=="log")
    family="poisson"
  if(fam$family=="gaussian" & fam$link=="identity")
    family="gaussian"
  if(fam$family=="Gamma" & fam$link=="log")
    family="gamma"
}

if ( is.character(family) ) {
    if (substr(family,1,3) == "gau") {
       familynum <- 0 # gaussian
       familyname <- "gaussian"
       linkfun <- 0 # not meaningful
    }
    else if (substr(family,1,3) == "poi") {
       familynum <- 1   #poisson
       familyname <- "poisson"
       linkfun <- 0 # not meaningful
    }
    else if (substr(family,1,3) == "neg") {
       familynum <- 2  #neg bin
       familyname <- "negative.binomial"
       linkfun <- 0 # not meaningful
    }
    else if (substr(family,1,3) == "bin") {
       familynum <- 3  #logistic/binomial
       familyname <- "binomial(link=logit)"
       link = "logit"  #by default
       linkfun <- 0    #logit link
    }
    else if (substr(family,1,3) == "clo") {
       familynum <- 3  #binomial(cloglog)
       familyname <- "binomial(link=cloglog)"
       link = "cloglog"
       linkfun <- 1
       K = 1 # response must be binary
    }
    else if (substr(family,1,3) == "gam") {
       familynum <- 4
       familyname<- "gamma"
       linkfun <- 0 # not meaningful, we only use the lo link at the moment
    }
    # FAMILY EDIT
    else stop (paste("'family'", family, "not recognised. See ?manyglm for currently available options."))
}
else stop("'family' not recognised. See ?manyglm for currently available options.")

#stop ("Current manyglm only supports the following link function for binary binomial regression: 'logit', 'cloglog'.")
ret.x <- x
ret.y <- y
ret.qr <- qr
cl <- match.call()
mf <- match.call(expand.dots = FALSE)
m  <- match(c("formula", "data", "subset", "na.action", "offset"), names(mf), 0)
mf <- mf[c(1, m)]

mf$drop.unused.levels <- FALSE
mf[[1]] <- as.name("model.frame")
mf <- eval(mf, parent.frame())    # Obtain the model.frame.
if(missing(data)) # Only coerce to model frame if not specified by the user.
  data <- mf

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

if (any(!is.wholenumber(Y)) & familynum != 4)
  warning(paste("Non-integer data are fitted to the", familyname, "model."))

# begin DW edits, 10/4/19
# if composition=TRUE, put in long format and send back to manyglm
if(composition==TRUE)
{
  # DW edits, 22/3/21, fixing issue #92
  # put in long format
  if(names(mf[1])%in%names(data)) #if response is in data frame
  {
    #DW edits, 1/2/22, so works when data is a list (like Code Box 14.6 of Eco-Stats text)
    whichResponse=which(names(data)==names(mf[1]))
    if(inherits(data,"data.frame")==FALSE)
      data2 = as.data.frame(data[-whichResponse])
    dat = data.frame(c(Y), data2[rep(1:N,p),])
    names(dat)[1] = names(mf)[1]  #match name to original object
    whichResponse=1
    # end DW edits
  }
  else
  {
      if(inherits(data,"data.frame")==FALSE)
        data = as.data.frame(data)
      dat = data.frame(c(Y), data[rep(1:N,p),])
      names(dat)[1] = names(mf)[1]  #match name to original object
      whichResponse=1
  }
  # make rows (row labels) and cols
  dat$rows = factor(rep(rownames(Y),p))
  dat$cols = factor(rep(colnames(Y),each=N))
  offset = rep(as.vector(model.offset(mf)),p) 
  #check no dollar signs in formula or database, these screw everything up
  names(dat)=gsub("$",".",names(dat),fixed=TRUE)
  
  # get formula for long format with composition
  if(length(mf)==1) #if no predictors, write formula with no cols interaction:
  {
    formLong=formula
    # now add cols, rows and cols:[prev formula]:
    formLong[3] = paste0(as.character(formula[3]),"+cols+rows")
    formLong=as.formula(paste0(formLong[2],formLong[1],formLong[3]))
  }
  else
  {
    if(formula[[3]]==".") #if ~., expand to variable names
      formula[[3]]=paste(names(mf[-whichResponse]),collapse="+")
    formLong=formula
    # now add cols, rows and cols:[prev formula]:
    formLong[3] = paste0("cols+",as.character(formula[3]),"+rows+cols:(",as.character(formula[3]),")")
    formLong[3] = gsub("$",".",formLong[3],fixed=TRUE)
    formLong=as.formula(paste0(formLong[2],formLong[1],formLong[3]))
  }
  z=manyglm(formLong, data=dat, block=dat$rows, composition=FALSE,
                 family=family, subset=subset, K=K, theta.method=theta.method,
                 model=model, x=x, y=y, qr=qr, cor.type=cor.type, 
                 shrink.param=shrink.param, tol=tol, maxiter=maxiter, 
                 maxiter2=maxiter2, show.coef=show.coef, show.fitted=show.fitted,
                 show.residuals=show.residuals, show.warning=show.warning, 
                 offset=offset,... )
    return(z)
  } #end David edits, 10/4/19
else{
  
  offset <- as.vector(model.offset(mf))
  
# HELP this doesn't make sense, was there previously a familynum == 4 that was meaningful
if ( familynum==3) {
    if(!is.null(labAbund) & all(substr(labAbund, 1,4)%in%c("succ", "fail")))
    {
       p <- p/2
       labAbund <- labAbund[substr(labAbund, 1,4) == "succ"]
       if(length(labAbund)!=p)
          stop( "for each variable of the response a column of successes and ",
                "a column of failures \nmust be provided ",
                "(marked by 'succ' and 'fail', see 'binstruc')")
    }
    if (all(is.wholenumber(Y)) & (length(Y[Y>1]>0)))
        warning("Count data are fitted to the binomial regression. Conisder a transformation first.")
    if ( (length(Y[Y<0])>0) | (length(Y[Y>1]>0)) )
    {
        if (familynum==3)
            warning("Data exceeds the range [0, 1].")
        if (familynum==4)
            stop("Data exceeds the allowed range [0, 1].")
    }
}

    ##################### BEGIN Estimation ###################
    # Obtain the Designmatrix.
    X <- model.matrix(mt, mf)

# Obtain regression parameters
if (familyname=="gaussian") {
    stop("Please use manylm to fit Gaussian")
#   z <- manylm(formula, data=data, subset=subset, na.action=na.action, model=model, x=x, y =y, qr=qr, cor.type=cor.type, shrink.param=shrink.param, tol=tol, ...)
#    z$family <- "gaussian"
#    z$formula <- formula
#    z$data <- data
#    class(z) <- c("manylm", "mlm")
#    return(z)
}
else {
    assign <- attr(X, "assign")
    tX <- t(X)
    dup <- duplicated(tX)
    assign <- assign[!dup]
    tX <- t(tX[!dup, ])  # remove duplicated col
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
    if (theta.method == "ML") methodnum <- 0
    else if (theta.method == "Chi2") methodnum <- 1
    else if (theta.method == "PHI") methodnum <- 2
    else if (theta.method == "MM") methodnum <- 3
    else stop("Check the param estimation method name.")

    if(family == 'gamma') {
        if (!(theta.method %in% c("ML", "MM"))) {
            theta.method <- "ML"
            methodnum <- 0
        }
    } else if(theta.method == "MM") stop("theta.method = 'MM' is only allowed for the gamma family.")

    if (show.warning==TRUE) warn <- 1
    else warn <- 0

    ######### call Glm Fit Rcpp #########
    modelParam <- list(tol=tol, regression=familynum, link=linkfun, estimation=methodnum, stablizer=FALSE, n=K, maxiter=maxiter, maxiter2=maxiter2, warning=warn)
    if(is.null(offset)) O <- matrix(0, nrow=N, ncol=p)
    else if (NCOL(offset)==1) O <- matrix(rep(offset), nrow=N, ncol=p)
    else O <- as.matrix(offset)
#    z <- .Call("RtoGlm", modelParam, Y, X, O, PACKAGE="mvabund")
    z <- RtoGlm(modelParam, Y, X, O)

    if(any(z$var.est==0)) #DW change, 30/10/14
    {
      z$var.estimator = pmax(z$var.est,1.e-6) #1.e-6 is used because it is elsewhere in this matrix. But shouldn't it be tol?
      z$residuals = (Y - z$fit) / sqrt(z$var.est)
    }

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
             warning("the calculated p-values might be unreliable as the number of cases is not much larger than the number of variables")
         shrink.param <- 1
    }

    # parameters that are needed for tests in anova.manylm and summary.manylm.
    dimnames(z$fitted.values) <- list(labObs, labAbund)
    dimnames(z$coefficients) <- list(colnames(X), labAbund)
    dimnames(z$var.coefficients) <- list(colnames(X), labAbund)
    dimnames(z$linear.predictor) <- list(labObs, labAbund)
    dimnames(z$residuals) <- list(labObs, labAbund)
    dimnames(z$PIT.residuals) <- list(labObs, labAbund)
    dimnames(z$sqrt.1_Hii) <- list(labObs, labAbund)
    dimnames(z$var.estimator) <- list(labObs, labAbund)
    dimnames(z$sqrt.weight) <- list(labObs, labAbund)
    names(z$theta) <- labAbund
    names(z$two.loglike) <- labAbund
    names(z$deviance) <- labAbund
    names(z$aic) <- labAbund
    names(z$iter) <- labAbund

    z$data <- data
    z$stderr.coefficients <- sqrt(z$var.coefficients)
    dimnames(z$stderr.coefficients) <- list(colnames(X), labAbund)
    z$phi <- 1/z$theta
    z$tol <- tol
    z$maxiter <- maxiter
    z$maxiter2 <- maxiter2
    z$prior.weight <- NULL
    z$AICsum <- sum(z$aic)
    z$family    <- familyname
    z$K         <- K
    z$theta.method<- theta.method
    z$cor.type  <- cor.type
    z$shrink.param  <- shrink.param
    z$call      <- cl
    z$terms     <- mt
    z$assign    <- assign
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
    z$offset <- O
    class(z) <- c("manyglm", "mglm")
    return(z)
  }
} # end composition!=TRUE
} #end function


is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)
{
    return(abs(x - round(x)) < tol)
}
