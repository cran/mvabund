predict.manylm <-
function (object, newdata, se.fit = FALSE, scale = NULL, df = Inf, 
    interval = c("none", "confidence",  "prediction"), level = 0.95,
    type = c("response", "terms"), terms = NULL, na.action = na.pass,
    pred.var = res.var/weights, weights = 1, ...) {


#    n <- NROW(object$residuals)
    p <- object$rank
    p1 <- seq_len(p)
    mf <- model.frame(object, data=object$data)
    nvar <- NCOL(model.response(mf))
    ynames <- colnames(as.matrix(model.response(mf)))

    tt <- terms(object)
    if (missing(newdata) || is.null(newdata)) {
        n <- NROW(object$residuals)
        mm <- X <- model.matrix(object, data=object$data)
        mmDone <- TRUE
        offset <- object$offset
    } else {
        n <- NROW(newdata) 
        Terms <- delete.response(tt)
#        m <- model.frame(Terms, newdata, na.action = na.action, 
#            xlev = object$xlevels)
#        if (!is.null(cl <- attr(Terms, "dataClasses"))) 
#            .checkMFClasses(cl, m)
#        X <- model.matrix(Terms, m, contrasts.arg = object$contrasts)
        X <- model.matrix(object)[1:n,]
        offset <- if (!is.null(off.num <- attr(tt, "offset"))) 
                       eval(attr(tt, "variables")[[off.num + 1]], newdata)
                  else if (!is.null(object$offset)) 
                       eval(object$call$offset, newdata)
        mmDone <- FALSE
    }
    
    isGlm <- inherits(object, "manyglm")
    
#    if(isGlm)    piv <- object$qr[[1]]$pivot[p1] else {
    piv <- object$qr$pivot[p1] #}

    if (p < ncol(X) && !(missing(newdata) || is.null(newdata))) 
        warning("prediction from a rank-deficient fit may be misleading")
    beta <- as.matrix(object$coefficients)
    
    predictor <- X[, piv, drop = FALSE] %*% beta[piv,, drop=FALSE]

    if (!is.null(offset)) 
        predictor <- predictor + offset
    interval <- match.arg(interval)
    if (interval == "prediction") {
        if (missing(newdata)) 
            warning("Predictions on current data refer to _future_ responses\n")
        if (missing(newdata) && missing(weights)) {
            w <- naresid(object$na.action, object$weights)
            if (!is.null(w)) {
                weights <- w
                warning("Assuming prediction variance inversely proportional to weights used for fitting\n")
            }
        }
        if (!missing(newdata) && missing(weights) && !is.null(object$weights) &&
            missing(pred.var)) 
            warning("Assuming constant prediction variance even though model fit is weighted\n")
            
        if (inherits(weights, "formula")) {
            # in contrast to predict.lm this is not allowed here
            stop("'weights' inherited from formula not yet implemented")
            
            if (length(weights) != 2L) 
                stop("'weights' as formula should be one-sided")
            d <- if(missing(newdata) || is.null(newdata))
                as.data.frame(mf)
            else newdata
            weights <- eval(weights[[2L]], d, environment(weights))
        }
        
    }
    type <- match.arg(type)
    if (se.fit || interval != "none") {

        res.var <- if (is.null(scale)) {           # nvar x nvar matrix
            r <- as.matrix(object$residuals)
            w <- object$weights
            if(!isGlm) {

              rss <- colSums(if (is.null(w)) r^2 else r^2 * w)
            } else {
                  rss <- colSums(if (is.null(w)) r^2 else r^2 * w)

            }
            df <- n - p
            rss/df

        } else scale^2

        if (type != "terms") {
            if (p > 0) {
                if(!isGlm){
                XRinv <- if (missing(newdata) && is.null(w))     # n x p matrix
                  qr.Q(object$qr)[, p1, drop = FALSE]
                else X[, piv] %*% qr.solve(qr.R(object$qr)[p1,
                  p1])
                  # only diag values are regarded, not really multivariate
                  ip <- XRinv^2 %*% matrix(rep(res.var, each=p), nrow=p,
                      ncol=nvar)
                } else {
                ip       <- matrix(nrow=n, ncol=nvar)
                for(ivar in 1:nvar)  {
                  XRinv <- if (missing(newdata) && is.null(w))   # n x p matrix
                  qr.Q(object$qr[[ivar]])[, p1, drop = FALSE]
                else X[, piv] %*% qr.solve(qr.R(object$qr[[ivar]])[p1,
                  p1])

                  ip[,ivar] <- drop(XRinv^2 %*% rep(res.var[ivar], p))
                }
                }
                
            } else ip <- matrix(0, nrow=n, ncol=nvar)
        }
    }
    ### begin type = "terms"
    if (type == "terms") {
        if (!mmDone) {
            mm <- model.matrix(object, data=object$data)
            mmDone <- TRUE
        }
        aa <- attr(mm, "assign")
        ll <- attr(tt, "term.labels")
        hasintercept <- attr(tt, "intercept") > 0L
        if (hasintercept) 
            ll <- c("(Intercept)", ll)
        aaa <- factor(aa, labels = ll)
        asgn <- split(order(aa), aaa)
        if (hasintercept) {
            asgn$"(Intercept)" <- NULL
            if (!mmDone) {
                mm <- model.matrix(object, data=object$data)
                mmDone <- TRUE
            }
            avx <- colMeans(mm)
            termsconst <- colSums(avx[piv] %*% beta[piv,])
        }
        nterms <- length(asgn)
        
        if (nterms > 0) {
            predictor <- array( dim= c(NROW(X), nvar, nterms))
            dimnames(predictor) <- list(rownames(X), ynames, names(asgn) )
            if (se.fit || interval != "none") {
                ip <- array(dim= c( NROW(X), nvar, nterms) )

                dimnames(ip) <- list(rownames(X), ynames, names(asgn))
                if(!isGlm){
                    Rinv <- qr.solve(qr.R(object$qr)[p1, p1])
                } else {
                    Rinv <- NULL
                    for(ivar in 1:nvar)
                        Rinv[[ivar]] <- qr.solve(qr.R(object$qr[[ivar]])[p1, p1])
                }
            }
            if (hasintercept) 
                X <- sweep(X, 2L, avx, check.margin = FALSE)
            unpiv <- rep.int(0L, NCOL(X))
            unpiv[piv] <- p1
            for (i in seq.int(1L, nterms, length.out = nterms)) {
                iipiv <- asgn[[i]]
                ii <- unpiv[iipiv]
                iipiv[ii == 0L] <- 0L
                predictor[,, i] <- if (any(iipiv > 0L))
                  X[, iipiv, drop = FALSE] %*% beta[iipiv,]
                else 0
                if (se.fit || interval != "none") {
                  if(!isGlm) {
                  ip[,, i] <- if (any(iipiv > 0L)) {
                    as.matrix(X[, iipiv, drop = FALSE] %*% Rinv[ii, 
                      , drop = FALSE])^2 %*%
                        matrix(rep(res.var, each=p), nrow=p, ncol=nvar)
                        
                  } else 0
                  } else {
                    for(ivar in 1:nvar) {
                    ip[,ivar, i] <- if (any(iipiv > 0L)) {
                    as.matrix(X[, iipiv, drop = FALSE] %*% Rinv[[ivar]][ii,
                      , drop = FALSE])^2 %*% rep(res.var[ivar], each=p)
                   } else 0
                  }
                  }
               }
            }
            if (!is.null(terms)) {
                predictor <- predictor[, , terms, drop = FALSE]
                if (se.fit) 
                  ip <- ip[, , terms, drop = FALSE]
            }
            
        } else {
            predictor <- ip <- array(0, dim= c(n, nvar, 0))
        }
        attr(predictor, "constant") <- if (hasintercept) 
            termsconst
        else rep(0, times= nvar)
    }
    #### end if (type == "terms")
    
    if (interval != "none") {
        tfrac <- qt((1 - level)/2, df)

        hwid <- tfrac * switch(interval, confidence = sqrt(ip),
            prediction = sqrt(ip + matrix(rep(pred.var, each=n),nrow=n,
                        ncol=nvar)))

        if (type != "terms") {
            pred <-  array(dim= c(n, nvar, 3))
            pred[,,1] <- predictor
            pred[,,2] <- predictor + hwid
            pred[,,3] <- predictor - hwid
            predictor <- pred
            dimnames(predictor)[[3]] <- c("fit", "lwr", "upr")
            
        } else {
            lwr <- predictor + hwid
            upr <- predictor - hwid
        }
    }
    
    if (se.fit || interval != "none") 
        se <- sqrt(ip)
    if (missing(newdata) && !is.null(na.act <- object$na.action)) {
        predictor <- napredict(na.act, predictor)
        if (se.fit) 
            se <- napredict(na.act, se)
    }
    if (type == "terms" && interval != "none") {
        if (missing(newdata) && !is.null(na.act)) {
            lwr <- napredict(na.act, lwr)
            upr <- napredict(na.act, upr)
        }
        list(fit = predictor, se.fit = se, lwr = lwr, upr = upr, df = df,
            residual.scale = sqrt(res.var))
    }
    else if (se.fit) 
        list(fit = predictor, se.fit = se, df = df, residual.scale =
         sqrt(res.var))
    else predictor
}


