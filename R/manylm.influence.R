################################################################################
# manylm.influence: a function that returns the diagonal(s) of the hat matrix  #
# (if the model is a manyglm, hat is a matrix and each column is the           #
# diag(hat matrix) * dispersion of the respective column of y,                 #
# the change in the estimated coefficients which results when the i-th case    #
# is dropped from the regression,                                              #
# a vector whose i-th element contains the estimate of the residual standard   #
# deviation obtained when the i-th case is dropped from the regression.        #
# and a vector of weighted (or for class glm rather deviance) residuals.       #
################################################################################

manylm.influence <- function (model, do.coef = TRUE) {

     # For the calculation of '.Fortran("lminfl",... ' : package "base" is safer
     # however, in R2.6.1 'Fortran symbol name "lminfl" is not in DLL for package "base" '
    if (R.version$major == "2" & R.version$minor == "6.1") {
      pack <- "stats"
    } else if (R.version$major == "2" & R.version$minor == "2.1") {
      pack <- "base"
    } else pack <- "base"      

    wt.res <- as.matrix(weighted.residuals(model))

    isManyglm <- inherits(model, "manyglm")

    # Avoid error for na.action = "na.exclude"
    e <- wt.res <- na.omit(wt.res)
    k <- as.integer(model$rank)
    # Calculate n, and NOTE: if y/x contained NA's the rows with NA's are not included
    n <- nrow(wt.res)
    if(any(attr(model$na,"class")=="omit") ) {  na.adjust <- length(model$na)
     } else  na.adjust <- 0
    # n.total <- n + na.adjust

    p <- ncol(wt.res)

    if (k == 0) {

        sigma <- sqrt(deviance(model)/df.residual(model))
        if(! isManyglm) {
        
            res <- list( hat = rep(0, n), coefficients = matrix(0, n, 0),
                sigma = matrix(rep(sigma, n), nrow=n, ncol=p, byrow=TRUE),
                wt.res = e )
                
        } else {
                
            hat <- matrix(0, nrow = n, ncol = p) 
            
            coefficients <- vector("list", p)
            
            if (do.coef) {
                for ( i in 1:p) coefficients[[i]] <- matrix(0, n, 0)
                names(coefficients ) <- colnames(coef(model))           
            }
            res <- list( hat = hat, coefficients = coefficients,
                sigma = matrix(rep(sigma, n), nrow=n, ncol=p, byrow=TRUE),
                wt.res = e )
        }
        
    } else {
        e[abs(e) < 100 * .Machine$double.eps * median(abs(e))] <- 0


        do.coef <- as.logical(do.coef)
        
        if(!isManyglm) {
 		     Qr <- model$qr

       	  n <- as.integer(nrow(Qr$qr))
        	k <- as.integer(Qr$rank)     
                
        if (NROW(e) != n) 
            stop("non-NA residual length does not match cases used in fitting")

        coefficients  <- vector("list", p)
        sigma         <-  matrix(nrow=n, ncol=p)
        
          for ( i in 1:p) {

            res <- try(.Fortran("lminfl", Qr$qr, n, n, k, as.integer(do.coef),
                Qr$qraux, wt.res = e[,i, drop=FALSE], hat = double(n),
                coefficients = if (do.coef) matrix(0, 
                    n, k) else double(0), sigma = double(n), tol = 10 * 
                    .Machine$double.eps, DUP = FALSE, PACKAGE = pack)[c("hat",  
                "coefficients", "sigma", "wt.res")] , silent=TRUE)
                
            if(inherits(res, "try-error"))  {
            # If the pack="base" and that didn't work, try package = "stats.
                res <- try(.Fortran("lminfl", Qr$qr, n, n, k, as.integer(do.coef),
                    Qr$qraux, wt.res = e[,i, drop=FALSE], hat = double(n),
                    coefficients = if (do.coef) matrix(0, 
                        n, k) else double(0), sigma = double(n), tol = 10 * 
                        .Machine$double.eps, DUP = FALSE, PACKAGE = "stats")[c("hat",
                    "coefficients", "sigma", "wt.res")] , silent=TRUE)

            if(inherits(res, "try-error")) stop(geterrmessage())
            
            }  
            coefficients[[i]] <- res$coefficients
            if (do.coef) names(coefficients ) <- colnames(coef(model))
            sigma[,i] <- res$sigma
            wt.res[,i] <- res$wt.res
         }
         res$coefficients   <- coefficients
         res$sigma          <- sigma
         res$wt.res         <- wt.res
                
        } else {

            # If isManyglm, model$qr is a list.
            Qr <- model$qr[[1]] 
        	  n <- as.integer(nrow(Qr$qr))
        	  k <- as.integer(Qr$rank)
                
        	if (NROW(e) != n) 
            stop("non-NA residual length does not match cases used in fitting")

            coefficients <- vector("list", p)
            hat <- matrix(0, nrow = n, ncol = p)
            sigma   <-  matrix(nrow=n, ncol=p)
            
            for ( i in 1:p) {

            Qr <- model$qr[[i]] 
            res <- try( .Fortran("lminfl", Qr$qr, n, n, k, as.integer(do.coef), 
                Qr$qraux, wt.res = e[,i, drop=FALSE], hat = double(n),
                coefficients = if (do.coef) matrix(0, 
                    n, k) else double(0), sigma = double(n), tol = 10 * 
                    .Machine$double.eps, DUP = FALSE, PACKAGE = pack)[c("hat",
                "coefficients", "sigma", "wt.res")], silent=TRUE)
            
            if(inherits(res, "try-error"))  {
            # If the pack="base" and that didn't work, try package = NULL.
                res <- try(.Fortran("lminfl", Qr$qr, n, n, k, as.integer(do.coef), 
                    Qr$qraux, wt.res = e[,i, drop=FALSE], hat = double(n),
                    coefficients = if (do.coef) matrix(0, 
                        n, k) else double(0), sigma = double(n), tol = 10 * 
                        .Machine$double.eps, DUP = FALSE, PACKAGE = "stats")[c("hat",
                    "coefficients", "sigma", "wt.res")] , silent=TRUE)
            
            } 
            
            hat[,i] <- res$hat          
            coefficients[[i]] <- res$coefficients   
            if (do.coef) names(coefficients ) <- colnames(coef(model))
            sigma[,i] <- res$sigma
            wt.res[,i] <- res$wt.res

            }
            # hat and coefficients must be calculated for all elemens of  model$qr.
            res$hat            <- hat
            res$coefficients   <- coefficients
            res$sigma          <- sigma
            res$wt.res         <- wt.res

        }

        if (!is.null(model$na.action)) {
            hat <- naresid(model$na.action, res$hat)
            hat[is.na(hat)] <- 0
            res$hat <- hat

            if (do.coef) {
                for ( i in 1:p) {
                    coefficientsi <- naresid(model$na.action, res$coefficients[[i]])
                    coefficientsi[is.na(coefficientsi)] <- 0
                    res$coefficients[[i]] <- coefficientsi
                }
            } 

            sigma <- naresid(model$na.action, res$sigma)
            sigma[is.na(sigma)] <- sqrt(deviance(model)/df.residual(model))
            res$sigma <- sigma
        }
    }
    res$wt.res <- naresid(model$na.action, res$wt.res)
    res$hat[res$hat > 1 - 10 * .Machine$double.eps] <- 1

    colnames(res$sigma) <-  colnames(res$wt.res)    

    if (!do.coef) {
            res$coefficients <- NULL
    } else {
          for ( i in 1:p) {
          rownames(res$coefficients[[i]]) <- rownames(res$wt.res) }
          if(ncol(res$coefficients[[1]])>0){
              coefs <- coef(model)
              nas <- rep(TRUE, times=NROW(coefs ) )
              for (i in 1:NROW(coefs ) )
                nas[i]<- ! any(is.na(coefs[i,]))
              for ( i in 1:p) {
              colnames(res$coefficients[[i]]) <- (rownames(coefs))[nas] }
          }
    }

    if(!isManyglm) {
        names(res$hat) <- rownames(res$sigma) <- rownames(res$wt.res) 
    } else {

        rownames(res$hat) <- rownames(res$sigma) <- rownames(res$wt.res) 
        colnames(res$hat) <- colnames(res$wt.res)
    }
    
    return(res)
}


