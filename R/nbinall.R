################################################################################
# NBINALL Log-likelihood function for a negative binomial, given MU.
# quadratic parameterization: Var(X)=E(X)+A.EX^2
# Can handle multivariate data
# for internal use
################################################################################

nbinall <- function( a, data, mu, lgamma= lgamma(data+1), tollevel=1.e-4 , 
	weights=rep.int(1, times = nRows))
{
	data 	<- as.matrix(data)
	lgamma	<- as.matrix(lgamma)
	nRows	<- nrow( data )
	nCols	<- ncol( data )
	
	# for mu == 0 the likelihood is asumed to be 0
  like    <- matrix( 0, nRows, nCols )
		
	if (length(mu)==1)  mu <- matrix(rep(mu,times=nRows*nCols),nrow=nRows, ncol=nCols)
		else 			mu <- as.matrix(mu)

  # avoid warning for lgamma([])
	isMuNon0 <- ( mu > tollevel )
 if(any(is.na(isMuNon0)))  print(mu)
	isPoi <- a == 0 
	
	isMuNon0.poi 		 <- isMuNon0
	isMuNon0.poi[, !isPoi] <- FALSE
	
	# find log-likelihood in Poisson case (excluding mu=0)
	if (any(isMuNon0.poi)) {
		like[isMuNon0.poi] <- data[isMuNon0.poi]*log( mu[isMuNon0.poi] ) -
      mu[isMuNon0.poi] - lgamma[isMuNon0.poi]
	}
	
	isMuNon0[, isPoi] <- FALSE
	
	# find LL in negbin case, where mu > tollevel
	if ( any(isMuNon0) )  {
	
		a 	<- matrix(rep(a, each = nRows), nrow=nRows, ncol= nCols)
		p   <- 1 / ( 1 + a*mu )	
		
		like[isMuNon0] <- lgamma( data[isMuNon0] + 1/a[isMuNon0] ) - lgamma(1/a[isMuNon0]) + 
				log(p[isMuNon0])/a[isMuNon0] + data[isMuNon0]*log(1-p[isMuNon0]) - lgamma[isMuNon0]
	}

	# corresponds to like <- -sum(w*like)
	# For multivariate data, To maximise likelihood, minimise -logL
	like <- - matrix(weights, nrow=1, ncol=nRows) %*% like 
	return(like)
}

