###################################################################
# loglik: a function that returns the Log-likelihood function 	  #
# for different families							                            #
###################################################################
# it is not checked whether the values of data, mu and disp are allowed
# if this is not ensured before-hand, nonsense-results or errors may result

loglik <- function( data, mu, family.char, disp, tollevel=1.e-4 ,
    weights=rep.int(1, times = nRows), omit=NULL, n=1, summed = TRUE) {

data    <- as.matrix(data)

if(!is.null(omit)) {
  data <- data[- omit,, drop=FALSE]
  weights <- weights[-omit]
}

nRows   <- nrow( data )
nCols   <- ncol( data )

if (length(mu)==1) {
  mu <- matrix(rep(mu,times=nRows*nCols),nrow=nRows, ncol=nCols)
} else {
    mu <- as.matrix(mu)
    if(nrow(mu)< nRows & (nrow(mu)+length(omit)== nrow(mu)))
          mu <- mu[- omit,, drop=FALSE]
}

if(family.char == "binomial") {            # correct binomial
    # n can be a vector!
    if(length(n ==1)) { n <- matrix(n, ncol=nCols, nrow= nRows)
    } else if (length(n)==nRows)  n <- matrix(rep(n, times=nCols), ncol=nCols,
                                        nrow= nRows)

    version2 <- TRUE
    if(version2){    # data and mu is percentage of successes
           data <- data*n
           mu   <- mu*n
           weights <- c(weights/n[,1])
    }   # else # n is no of tries, data is number of successes

    if(any(data==0))   stop("data must be >0")
    if(any(data>n))   stop("data must <n")

    like <- lgamma(n+1) - lgamma(data + 1) -lgamma(-data + n + 1) +
            data*log(mu/n) +
            # ( mu-data)*log(1- mu/data)  # version in Hardin & Hilbe
            ( n-data)*log(1- mu/n)    # seems to be the right one

    if(!summed) {return( weights*(like) ) }  else {
    like <-  c( matrix(weights, ncol=nRows, nrow=1) %*% (like) )   # disp = 1 for binomial
    return(c(like)) }

} else if(family.char  == "gaussian"){

    like <- (data - mu)^2
    like <- like /(matrix(rep(disp, each=nRows),nRows,nCols))
    like <- like + matrix( rep(log(2*pi*disp), each=nRows),nRows,nCols )

    if(!summed) {
       return( weights*(- 1/2)*like )
    } else {
        like <- -1/2 * like
        like <- c( matrix(weights, ncol=nRows, nrow=1) %*% (like) )
        return(c(like))
    }
} else if(family.char  == "Gamma"){
    
    like <- data/mu - log(matrix(rep(disp, each=nRows), nrow= nRows, ncol= nCols )/mu)
 
    if(!summed) {
       like <- like - matrix(rep((disp - 1)/disp, each=nRows),nRows,nCols)*log(data)+
        matrix(rep(lgamma(disp)/disp, each=nRows),nRows,nCols)
      return( weights*matrix(rep(-1/disp, each=nRows),nRows,nCols)*like )
    } else {
      like <- like - matrix(rep((disp - 1)/disp, each=nRows),nRows,nCols)*log(data) +
        matrix(rep(lgamma(disp)/disp, each=nRows),nRows,nCols)
      like <- -1/disp * c( matrix(weights, ncol=nRows, nrow=1) %*% like )
      return(c(like))
    }
 
} else if(family.char  == "inverse.gaussian") {
     like1 <- log(data^3 * matrix(rep(disp, each = nRows), nrow=nRows, ncol=nCols))
     like2 <- ( (data - mu)^2/(data*mu^2 ) ) /
              matrix(rep(disp, each=nRows),nRows,nCols)
     like <- like1 + like2 + matrix(rep(nRows * log(2*pi), each=nRows),nRows,nCols)

    if(!summed) {
      like <- -1/2 * weights *like
      return(like)
    } else {
      like <- c( -1/2 * matrix(weights, ncol=nRows, nrow=1) %*% like )
      return(c(like))
    }

} else if(family.char  == "poisson") {

    like <- - mu + data*log(mu) - lgamma(data+1)
    if(!summed) {return( weights*(like) ) }  else {
    like <-  c( matrix(weights, ncol=nRows, nrow=1) %*% (like)  )
    # disp = 1 for poisson
    return(like) }

}  else if(family.char  == "quasipoisson") {
    # disp must be available 'a priori' to calc ll
    like     <-  matrix( 0, nRows, nCols )
    isMuNon0 <- ( mu > tollevel )
    isPoi <- disp == 1

	 isMuNon0.poi 		 <- isMuNon0
	 isMuNon0.poi[, !isPoi] <- FALSE

  	# find log-likelihood in Poisson case (excluding mu=0)
	 if(any(isMuNon0.poi)) {
          like[isMuNon0.poi] <- data[isMuNon0.poi]*log( mu[isMuNon0.poi] ) -
            mu[isMuNon0.poi] - lgamma(data+1)[isMuNon0.poi]
	 }

	 isMuNon0[, isPoi] <- FALSE		# case isPoi already done, now ignore these

    if ( any( isMuNon0 ) )    {
		a 		<- matrix(rep(disp, each = nRows), nrow=nRows, ncol= nCols)
		p		<- 1 / a
		r		<- mu / ( a - 1 )
   like[isMuNon0] <- lgamma( data[isMuNon0] + r[isMuNon0]) -
      lgamma( r[isMuNon0] ) + r[isMuNon0]*log(p[isMuNon0]) +
      data[isMuNon0]*log(1-p[isMuNon0])- lgamma(data+1)[isMuNon0]
	}

  if(!summed) {return( weights*(like) ) }  else {
	like <- matrix(weights, nrow=1, ncol=nRows) %*% like
  return(c(like)) }

} else if(family.char  == "negative.binomial" | family.char  ==
    "Negative Binomial") {  # disp must be available 'a priori'  to calc ll
    like     <-  matrix( 0, nRows, nCols )
    isMuNon0 <- ( mu > tollevel )
    isPoi <- disp == 0

  	isMuNon0.poi 		 <- isMuNon0
  	isMuNon0.poi[, !isPoi] <- FALSE

  	# find log-likelihood in Poisson case (excluding mu=0)
  	if (any(isMuNon0.poi)) {
  		like[isMuNon0.poi] <- data[isMuNon0.poi]*log( mu[isMuNon0.poi] ) -
        mu[isMuNon0.poi] - lgamma(data+1)[isMuNon0.poi]
  	}

  	isMuNon0[, isPoi] <- FALSE		# case isPoi already done, now ignore these

    	# find LL in negbin case, where mu > tollevel
	  if ( any(isMuNon0) )  {

		a 	<- matrix(rep(disp, each = nRows), nrow=nRows, ncol= nCols)
		p   <- 1 / ( 1 + a*mu )

		like[isMuNon0] <- lgamma( data[isMuNon0] + 1/a[isMuNon0] ) -
      lgamma(1/a[isMuNon0]) + log(p[isMuNon0])/a[isMuNon0] +
      data[isMuNon0]*log(1-p[isMuNon0]) - lgamma(data+1)[isMuNon0]
	}

  if(!summed) {return( weights*(like) ) }  else {
	like <- matrix(weights, nrow=1, ncol=nRows) %*% like
  return(c(like))  }

} else if(family.char == "quasibinomial") {
  # at the moment this is the same as binomial, corrections necessary

    if(length(n ==1)) { n <- matrix(n, ncol=nCols, nrow= nRows)
    } else if (length(n)==nRows)  n <- matrix(rep(n, times=nCols),
      ncol=nCols, nrow= nRows)

    version2 <- TRUE
    if(version2){    # data and mu is percentage of successes
           data <- data*n
           mu   <- mu*n
    }   # else # n is no of tries, data is number of successes

    like <- lgamma(n+1) - lgamma(data + 1) -lgamma(-data + n + 1) +
            log(mu/n) +
            (data-1)*log(disp*data + mu/n) +( n-data)*log(1 - mu/n- disp*data )

    if(!summed) {return( weights*(like) ) }  else {
    like <-  c( matrix(weights, ncol=nRows, nrow=1) %*% (like) )
    return(c(like))  }

} else if(family.char  == "quasi") {
    warning("likelihood not yet implemented for family quasi")
    if(!summed) {return( matrix(nrow=NROW(data), ncol=NCOL(data)) ) }  else {
    return(rep(NA, times=NCOL(data)))  }
  }

} 


