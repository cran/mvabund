genRNDfiles <- function(object, resample=object$resample, nBoot=object$nBoot, filename=NULL)
{ 
   # generate random sequence and save in the corresponding file 
   nRows <- NROW(object$y)   # number of site observations 
   paths <- .find.package("mvabund2")
   if (resample == "score") {
      if ( filename==NULL )
         filename <- file.path(paths, "data", "scores.dat")
      perm  <- matrix(rnorm(nBoot*nRows),nrow=nBoot , ncol=nRows)
   }
   else if (resample == "perm.resid") {
      if (filename==NULL)
         filename <- file.path(paths, "data", "permID.dat")
      perm <- matrix( 0, nBoot, nRows )
      # random number generator from uniform distribution
      runifs  <- matrix( runif((nBoot-1)*nRows),nrow=(nBoot-1),ncol=nRows)		
      for (i in 1:nBoot-1)
          perm[i,] <- order( runifs[i,] ) 
      perm[nBoot, ] <- c(1:nRows)
   }
   else if (resample == "case" | resample=="residual") {
      if (filename==NULL)
         filename <- file.path(paths, "data", "bootID.dat")
      perm     <- ceiling( nRows*runif( nBoot*nRows ) )
      dim(perm) <- c( nBoot, nRows )
   }
   else {   # for test purpose
      if (filename==NULL)
          filename <- file.path(paths, "data", "tmp.dat")
      perm  <- matrix(rnorm(10*nRows),nrow=10 , ncol=nRows)
   }
      
   write.table(perm, filename, col.names=F, row.names=F)

}
