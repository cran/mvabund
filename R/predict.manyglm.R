
predict.manyglm <- 
function (object, newdata = NULL, type = c("link", "response", 
    "terms"), se.fit = FALSE, dispersion = object$phi, terms = NULL, 
    na.action = na.pass, ...) {
    
    ny      <- NCOL(object$fitted.values)
    nobs    <- NROW(object$fitted.values)
    type    <- match.arg(type)
    na.act  <- object$na.action
    object$na.action <- NULL
#    newdata <- data.frame(newdata)

    if (!se.fit) {
        if (missing(newdata)) {
            pred <- switch(type, link = object$linear.predictor, 
                response = object$fitted.values, terms = predict.manylm(object,                se.fit = se.fit, scale = 1, type = "terms", terms = terms))
            if (!is.null(na.act)) 
                pred <- napredict(na.act, pred)
        }
        else {
            pred <- predict.manylm(object, newdata=newdata, se.fit, scale = 1, 
                type = ifelse(type == "link", "response", type), 
                terms = terms, na.action = na.action)
            switch(type, response = {
                pred <- matrix(family(object)$linkinv(pred),
                nrow = nobs , ncol= ny ) },
            link = , terms = )
        }
        
    } else {
#        if (inherits(object, "survreg")) 
##            dispersion <- rep(1, times=ncol(object$y))
#             dispersion <- dispersion
#        if (is.null(dispersion) || dispersion == 0) 
#            dispersion <- object$phi
#        if (is.null(dispersion) || dispersion == 0)
#            dispersion <- summary(object,dispersion = dispersion, nBoot=2)$phi
            # test values not interesting here, therefore set nBoot small for
            # faster calculation
        residual.scale <- as.vector(sqrt(dispersion))
        pred <- predict.manylm(object, newdata, se.fit, scale = residual.scale, 
            type = ifelse(type == "link", "response", type), 
            terms = terms, na.action = na.action)
        fit <- pred$fit
        se.fit <- pred$se.fit

        switch(type, response = {
            mu.eta.val <- matrix(family(object)$mu.eta(fit), nrow = nobs ,
                          ncol= ny )
            se.fit <- se.fit * abs(mu.eta.val)
            fit <- matrix(family(object)$linkinv(fit), nrow = nobs , ncol= ny )
                        
        }, link = , terms = )
        if (missing(newdata) && !is.null(na.act)) {
            fit <- napredict(na.act, fit)   # test with missing values in object
            se.fit <- napredict(na.act, se.fit)
        }
        pred <- list(fit = fit, se.fit = se.fit, residual.scale = residual.scale)
    }
    pred
}

