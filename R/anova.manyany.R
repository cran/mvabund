anova.manyany = function(object, ..., nBoot=99, p.uni="none", block = object1$block, bootID=NULL, replace=TRUE)
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
           warning("current version only compares two manyany objects")
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
    if(grepl("egative",object1$family[[i.var]]$family) || object1$family[[i.var]]$family == "negbinomial")
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
    if(object1$family[[i.var]]$family=="ordinal")
      qfn[i.var] = "qordinal"
  }

  if(is.null(bootID)==FALSE)
  {
    bootID = as.matrix(bootID)
    if(dim(bootID)[2]!=n.rows)
      stop("Number of rows of bootID must match number of rows in data")
    nBoot = dim(bootID)[1] #overwriting nBoot with value implied by user-entered ID matrix
    block = NULL #overwriting previous value for block
    print("User-entered bootID matrix will be used to generate bootstrap samples")
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
#  ft.1i=eval(object1$call) #this call seems unnecessary
#  ft.2i=eval(object2$call) #this call seems unnecessary
  statj = 2 * ( logLik(object2)-logLik(object1) )
  stat = sum(statj)
  
  #initialise parameters for bootstrapping
  yMat = matrix(NA,n.rows,n.vars)
  if(object1$family[[1]]$family=="ordinal")
    yMat=data.frame(yMat)
  statj.i = matrix(NA,n.vars,nBoot)
  if(n.vars>1)
    dimnames(statj.i)[[1]] = dimnames(object1$residuals)[[2]]
  stat.i=rep(NA,nBoot)
  if(is.null(bootID))
    boot.Resamp = rep(NA,n.rows)
  object1$call$get.what="none" #to avoid wasting time computing residuals etc when resampling
  object2$call$get.what="none" #to avoid wasting time computing residuals etc when resampling
    
  #now do the bootstrap
  for(iBoot in 1:(nBoot))
  {
    if(is.null(bootID)==FALSE)
      boot.Resamp = bootID[iBoot,]
    else      # generate resampled residuals
    {
      if(is.null(block))
        boot.Resamp = sample(1:n.rows,replace=replace)
      else
        boot.Resamp[unlistIDs] = unlist(blockIDs[sample(n.levels,replace=replace)]) #unlistIDs is needed to make sure each unlisted blockID ends up in the right place
    }
    resid.i = as.matrix(object1$residuals[boot.Resamp,])
    # simulate data to get resampled yMat

    for(i.var in 1:n.vars)
    {
      qparams = object1$params[[i.var]]
      qparams[[1]]=resid.i[,i.var]
      names(qparams)[1]="p"
      yMat[,i.var] = do.call(qfn[i.var], qparams)
    }
    #save resampled yMat as whatever the original yMat was called in workspace - but without zerotons
    if(object1$family[[1]]$family=="ordinal")
      is.zeroton = apply(yMat,2,function(x) length(table(x)))==1
    else
      is.zeroton = apply(yMat,2,sum)==0
    assign(as.character(object1$call[[3]]),yMat[,is.zeroton==FALSE]) 
    assign(as.character(object2$call[[3]]),yMat[,is.zeroton==FALSE]) 
    #re-fit manyany functions and calculate test stats using the resampled yMat:
    if(sum(is.zeroton==FALSE)>0)
    {
      ft.1i=eval(object1$call)
      ft.2i=eval(object2$call)
      statj.i[is.zeroton==FALSE,iBoot]=2 * ( logLik(ft.2i)-logLik(ft.1i) )
      stat.i[iBoot] = sum(statj.i[,iBoot], na.rm=TRUE)      
    }
    else
      stat.i[iBoot] = 0
  }
  p = ( 1 + sum(stat.i>stat-1.e-8) ) / (nBoot + 1)
  pj = ( 1 + apply(statj.i>statj-1.e-8,1,sum) ) / ( nBoot + 1)

  class(stat.i) = "numeric"
  if(p.uni=="unadjusted")
    result = list(stat=stat,p=p,uni.test=statj,uni.p=pj,stat.i=stat.i,statj.i=statj.i,p.uni=p.uni,nBoot=nBoot) 
  if(p.uni=="none")
    result = list(stat=stat,p=p,stat.i=stat.i,p.uni=p.uni,nBoot=nBoot) 

  class(result) = "anova.manyany"
  return(result)  
}

print.anova.manyany=function(x, ...)
{
  #get overall results in a table
  table=matrix(c(x$stat,x$p),1,2)
  dimnames(table)[[2]]=c("LR","Pr(>LR)")
  dimnames(table)[[1]]=c("sum-of-LR")

  allargs <- match.call(expand.dots = FALSE)
  dots <- allargs$...
  s.legend = TRUE
  if(length(dots)>1)
  {
    if("signif.legend" %in% dots)
      s.legend = signif.legend
  }
  if(x$p.uni=="none")
    signif.legend = s.legend
  else
    signif.legend = FALSE
  
  #print overall results
  cat("\n")
  printCoefmat(table, tst.ind=1, P.values=TRUE, has.Pvalue=TRUE, signif.legend=signif.legend, eps.Pvalue=1/(x$nBoot+1-1.e-8),...)
  cat("\n")
  #print univariate results in a table, if required
  if(x$p.uni!="none")
  {
    signif.legend = s.legend
    tablej=cbind(x$uni.test,x$uni.p)
    dimnames(tablej)[[2]]=c("LR","P(>LR)")
    printCoefmat(tablej, tst.ind=1, P.values=TRUE, has.Pvalue=TRUE, signif.legend=signif.legend, eps.Pvalue=1/(x$nBoot+1-1.e-8), ...)
  }
}

