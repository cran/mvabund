anova.manyany = function(object, ..., nBoot=99, p.uni="none", block = object1$block, replace=TRUE)
{
# analysis of variance comparing object1 (null) to object2 (alternative)
# uses the PIT-trap
# unadjusted univariate P-values only at this stage
# block is a factor specifying the sampling level to be resampled. For example, if data have multiple
#  rows of records for each site, e.g. multi-species data with entries for different species on different
#  rows, you can use your site ID variable as the block argument to resample sites only, for valid
#  cross-site inferences despite within-site species correlation. Well, valid assuming sites are independent.
#  You could do similarly for a repeated measures design to make inferences robust to temporal autocorrelation.
#  Note that block needs to be balanced to use a residual resampling approach, e.g. equal number of species
#  entries for each site (i.e. include rows for zero abundances too). Default is resampling rows (if
#  composition=TRUE in the manyany command, this means resampling rows of data as originally sent to manyany).
# replace=F for PIT-permutation, replace=T for PIT-trap
# 
# glmm fans - Sorry, we can't verify the method for random effects models yet
#  Probably needs resampling of random effects as well as observed data.
#  At best is approx valid for mixed models only for conditional inference (conditional on observed random effects)
#
# EXAMPLES
#
# ft=manyany("glm",abund,data=X,y~1,family="poisson")
# ftSoil=manyany("glm",abund,data=X,y~soil.dry,family="poisson")
# an=anova(ft,ftSoil,p.uni="unadjusted")
  
  object1 = object 
  # get object 2
    dots <- list(...)
    ndots <- length(dots)
    fndObj2 <- FALSE
    if (ndots==0) {
       stop("missing a second manyany object")
    }else {
       if (ndots>1)
           warning("crrently version only compares two manyany objects")
       for (i in 1:ndots) {
           if (any(class(dots[[i]])=="manyany")){
              object2 <- dots[[i]]
              fndObj2 <- TRUE
              break
           }
       }
       if (!fndObj2) stop("cannot find object 2")
    }

  if(any(names(object1$call)=="composition"))
  {
    if(object1$call$composition==TRUE) #recode so that it fits compositional models as univariate, to save time and fuss/bother.
    {
      object1$call$formula = object1$formula
      object2$call$formula = object2$formula
      object1$call$data = object1$model
      object2$call$data = object2$model
      object1$residuals = as.matrix(c(object1$residuals))
      object1$call$composition=FALSE
      object2$call$composition=FALSE
      assign(as.character(object1$call[[3]]),object1$model$y) 
      assign(as.character(object2$call[[3]]),object2$model$y) 
    }
  }
  n.rows = dim(object1$resid)[1]
  n.vars = dim(object1$resid)[2]
  
  qfn = rep(NA,n.vars)
  for(i.var in 1:n.vars)
  {
    if(object1$family[[i.var]]$family=="negbinomial")
      qfn[i.var] = "qnbinom"
    if(object1$family[[i.var]]$family=="poisson")
      qfn[i.var] = "qpois"
    if(object1$family[[i.var]]$family=="binomial")
    {
      qfn[i.var] = "qbinom"
      warning("The binomial option of manyany currently assumes you have presence/absence data")
    } 
    if(object1$family[[i.var]]$family=="gaussian")
      qfn[i.var] = "qnorm"  
    if(object1$family[[i.var]]$family=="Tweedie")
      qfn[i.var] = "qtweedie"
  }

  if(is.null(block)==FALSE)
  {
    tb=table(block)
    n.levels = length(tb)
    if(any(tb!=n.rows/n.levels))
    {   
      print(tb) 
      stop("Sorry, block needs to be a balanced factor - same number of rows for each level")
    }
    else
    {
      blockIDs = vector("list",n.levels)
      for(i.level in 1:n.levels)
        blockIDs[[i.level]] = which(block==names(tb)[i.level])
      unlistIDs = unlist(blockIDs) #needed to match each resampled observation with its correct location
    }
  }
  #get observed test stat
  ft.1i=eval(object1$call)
  ft.2i=eval(object2$call)
  statj = 2 * ( logLik(ft.2i)-logLik(ft.1i) )
  stat = sum(statj)
  
  #initialise parameters for bootstrapping
  yMat = matrix(NA,n.rows,n.vars)
  statj.i = matrix(NA,n.vars,nBoot)
  dimnames(statj.i)[[1]] = dimnames(object1$residuals)[[2]]
  stat.i=rep(NA,nBoot)
  boot.Resamp = rep(NA,n.rows)
  for(iBoot in 1:(nBoot))
  {
    # generate resampled residuals
    if(is.null(block))
      boot.Resamp = sample(1:n.rows,replace=replace)
    else
      boot.Resamp[unlistIDs] = unlist(blockIDs[sample(n.levels,replace=replace)]) #unlistIDs is needed to make sure each unlisted blockID ends up in the right place
      resid.i = as.matrix(object1$residuals[boot.Resamp,])
    # simulate data to get resampled yMat
    
    for(i.var in 1:n.vars)
    {
      qparams = object1$params[[i.var]]
      qparams[[1]]=resid.i[,i.var]
      names(qparams)[1]="p"
      yMat[,i.var] = do.call(qfn[i.var], qparams)
    }
    #save resampled yMat as whatever the original yMat was called in workspace
    assign(as.character(object1$call[[3]]),yMat) 
    assign(as.character(object2$call[[3]]),yMat) 

    #re-fit manyany functions and calculate test stats using the resampled yMat:
    ft.1i=eval(object1$call)
    ft.2i=eval(object2$call)
    statj.i[,iBoot]=2 * ( logLik(ft.2i)-logLik(ft.1i) )
    stat.i[iBoot] = sum(statj.i[,iBoot])
  }
  p = ( 1 + sum(stat.i>stat-1.e-8) ) / (nBoot + 1)
  pj = ( 1 + apply(statj.i>statj-1.e-8,1,sum) ) / ( nBoot + 1)

  class(stat.i) = "numeric"
  if(p.uni=="unadjusted")
    result = list(stat=stat,p=p,uni.test=statj,uni.p=pj,stat.i=stat.i,statj.i=statj.i,p.uni=p.uni) 
  if(p.uni=="none")
    result = list(stat=stat,p=p,stat.i=stat.i,p.uni=p.uni) 
  
  class(result) = "anova.manyany"
  return(result)  
}

print.anova.manyany=function(object)
{
  #get overall results in a table
  table=matrix(c(object$stat,object$p),1,2)
  dimnames(table)[[2]]=c("LR","P(>LR)")
  dimnames(table)[[1]]=c("sum-of-LR")
  
  #print overall results
  cat("\n")
  printCoefmat(table,has.Pvalue=T,signif.legend=F)
  cat("\n")
  
  #print univariate results in a table, if required
  if(object$p.uni!="none")
  {
  tablej=cbind(object$uni.test,object$uni.p)
  dimnames(tablej)[[2]]=c("LR","P(>LR)")
  printCoefmat(tablej,has.Pvalue=T)
  }
}

