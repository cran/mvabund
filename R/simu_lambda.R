################################################################################
##  simulation test: estimate shrink parameter on the generated data set simuY #
##  imported from "/c/z3069443/home/workbench/C/mvabund2/simuYdata/simuY_i.dat #
##  the esimated lambda for each simulated Y file is output to a single file:
##  "/c/z3069443/home/workbench/C/mvabund2/data/lambda.dat"
################################################################################

## import simuY

simu_lambda <- function( simuid, sigmaid ) {

    paths <-  "/c/z3069443/home/workbench/C/mvabund2"
    simuname <- paste("simu", simuid, sep="")
    sigmaname <- paste("sigma", sigmaid, sep="")
    outfile <- file.path(paths, "SimuYdata_log", simuname, sigmaname, "lambda_R.dat")
    # rm outfile
    sum <- 0
    for( i in 0:999 ) {
       fname <- paste("simuY_", i, ".dat", sep="")
       filename <- file.path(paths, "SimuYdata_log", simuname, sigmaname, fname)  
       simuY <- matrix(scan(filename, quiet=TRUE), nrow=12, ncol=28) 
       simuY <- t(simuY)
       shrink.param <- ridgeParamEst(dat=simuY, X=model.matrix(spiddat~spidx), only.ridge=TRUE)$ridgeParameter
       sum <- sum+shrink.param
#       cat(shrink.param, sep="\n")
       cat(shrink.param, file = outfile, sep = "\n", append = TRUE)
   }
   return (sum/1000)
}
