manyany = function(fn, yMat, formula, data, family="negative.binomial", composition = FALSE, var.power=NA, ...)
{
  #MANYANY - applies a function of your choice to each column of YMAT and computes logLik by taxon.
  # FN is a character vector giving the name of the function to be applied to each taxon.  e.g. "glm"
  # YMAT is a matrix contatining the response variable for each taxon.
  # FORMULA is the formula to use in the call to FM.
  # a FAMILY argument needs to be specified. It can be a list with different families for different variables (with length matching the number of columns in YMAT)
  # COMPOSITION is a logical switch indicating whether or not to do a compositional analysis (i.e., include a
  #  row effect in the model to account for changes in total abundance across samples).
  # VAR.POWER is needed for tweedie distributions - the power parameter needs to be specified here as well as in the family argument.
  # ... any further arguments required by FN.
  #
  # Examples:
  # require(mvabund)
  # data(spider)
  # abund=spider$abund
  # X=data.frame(spider$x)
  # 
  # a manygam:
  # library(mgcv)
  # ft=manyany("gam",abund,y~s(soil.dry),data=X,family="poisson")
  # 
  # a manyglmm:
  # library(lme4)
  # gr = rep(1:2,each=14)
  # ft=manyany("lmer",abund,y~soil.dry+1|gr,data=X,family="poisson")
  #
  ## A manyglm:
  # ft=manyany("glm",abund,y~soil.dry,data=X,family="poisson")
  ## note this gives the same answer as:  
  # ft=manyglm(mvabund(abund)~X$soil,family="poisson")
  
  yMat = as.matrix(yMat)
  data = as.data.frame(data)
  yNames = dimnames(yMat)
  n.rows = dim(yMat)[1]
  n.vars = dim(yMat)[2]
  if(is.null(dimnames(yMat)[[1]]))
    dimnames(yMat)[[1]] = 1:n.rows
  
  call=match.call()
  
  if(composition==FALSE)
  {
    formula = formula(paste("y~",formula[3],sep=""))
    block = NULL
  }
  else
  {
    yVec    = as.vector(yMat)
    ref     = factor(rep(1:n.rows,n.vars))
    spp     = factor(rep(1:n.vars,each=n.rows))
    data    = data.frame(ref,spp,data[ref,])
    formula = formula(paste("y~ref+spp+spp:(",formula[3],")",sep=""))
    n.rows.orig = n.rows #save for later
    n.vars.orig = n.vars #save for later
    n.rows  = length(yVec)
    n.vars  = 1
    names(yVec) = paste( dimnames(yMat)[[2]][spp], ".", dimnames(yMat)[[1]][ref], sep="")
    yMat    = as.matrix(yVec)
    if(class(family)!="family" & length(family)>1)
      stop("when using composition=TRUE, family argument must have length one.")
    block = ref #to make sure resampling is by row of original data, not of vectorised data.
  }

  #If family is specified once, turn it into a list of n.vars family objects
  if(class(family)=="family" || length(family)==1)
  {
    fam = family #temporary store to slot into a big list
    family = vector("list",n.vars)
    for(i.var in 1:n.vars)
      family[[i.var]] = fam
  }
  if(length(family)!=n.vars)
      stop("family argument has length more than one but not equal to the number of columns in yMat (!?)")

  if(length(var.power)==1)
    var.power=rep(var.power,n.vars)
    
  fam = family
  # now ensure each family is a proper family function
  for(i.var in 1:n.vars)
  {
    if (is.character(family[[i.var]])) 
    {
      if (family[[i.var]] == "negbinomial" || family[[i.var]]=="negative.binomial")
      {
        fam[[i.var]] = negative.binomial(10^6)
        fam[[i.var]]$family = family[[i.var]]
      }
      else
      {
        fam.fn = get(fam[[i.var]], mode = "function", envir = parent.frame())        
        fam[[i.var]] = fam.fn()
      }  
    }
    if(fam[[i.var]]$family=="binomial")
      warning("The binomial option of manyany currently assumes you have binary (presence/absence) response")
  }
  
  manyfit = list()
  fits = matrix(NA,n.rows,n.vars)
  etas = matrix(NA,n.rows,n.vars)
  params = list()
  logL = rep(NA,n.vars)
  for(i.var in 1:n.vars)
  {
    data$y = yMat[,i.var]
    manyfit[[i.var]] = do.call(fn, list(formula=formula, family=family[[i.var]], data=data, ...)) #note use of family argument as originally specified
    fits[,i.var] = fitted(manyfit[[i.var]])
    if(fn=="lmer")
      etas[,i.var] = manyfit[[i.var]]@eta
    else
      etas[,i.var] = predict(manyfit[[i.var]])
    logL[i.var] = logLik(manyfit[[i.var]])
    if(is.na(logL[i.var]))
       logL[i.var] = -0.5*deviance(manyfit[[i.var]]) #just in case logL function is undefined, e.g. tweedie 
    if(fam[[i.var]]$family == "negbinomial" || fam[[i.var]]$family=="negative.binomial")
    {
      if(any(names(manyfit[[i.var]])=="theta"))
        theta=manyfit[[i.var]]$theta
      else
        theta=1/manyfit[[i.var]]$phi
      params[[i.var]] = list(q=yMat[,i.var],mu=fits[,i.var],size=theta)
    }
    if(fam[[i.var]]$family=="poisson")
      params[[i.var]] = list(q=yMat[,i.var],lambda=fits[,i.var])
    if(fam[[i.var]]$family=="binomial")
      params[[i.var]] = list(q=yMat[,i.var],prob=fits[,i.var],size=1)
    if(fam[[i.var]]$family=="gaussian")
    {
      s.ft=summary(manyfit[[i.var]])
      if(any(names(s.ft)=="sigma"))
        sd=s.ft$sigma
      else
        sd=s.ft$scale
      params[[i.var]] = list(q=yMat[,i.var],mean=fits[,i.var],sd=sd)
    }
    if(fam[[i.var]]$family=="Tweedie")
    {
      params[[i.var]] = list(q=yMat[,i.var], power=var.power[i.var], mu=fits[,i.var], phi=summary(manyfit[[i.var]])$disp)
    }
  } 
  names(params) = dimnames(yMat)[[2]]
  attributes(logL)$df = attributes(logLik(manyfit[[i.var]]))$df
  attributes(logL)$nobs = n.rows
  names(logL) = dimnames(yMat)[[2]]
  class(logL) = "logLik"
  
  if(composition==TRUE) #reshape to original data size if required
  {
    fits   = matrix(fits,n.rows.orig,n.vars.orig)
    etas   = matrix(etas,n.rows.orig,n.vars.orig)
  }    
  resids = residuals.manyany(list(params=params, family=fam, composition=composition, fitted.values=fits))
  dimnames(resids) = yNames
  dimnames(fits)   = yNames
  dimnames(etas)   = yNames
  object=list(logL=logL,fitted.values=fits,residuals=resids,linear.predictor=etas,family=fam,call=call,params=params,model=model.frame(manyfit[[i.var]]), terms = terms(manyfit[[i.var]]), formula=formula, block=block, composition=composition)
  class(object)=c("manyany", class(manyfit[[i.var]]) )
  
  return(object)
}


print.manyany <- function(object, digits = max(3L, getOption("digits") - 3L),...)
{
  n.vars = dim(object$fitted)[2]
  cat("\nCall:  ", paste(deparse(object$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("Number of rows:\t   ", dim(object$fitted)[1], "\t Number of columns:\t   ", n.vars)
  cat("\n")
  cat("Number of parameters in model:\t   ", n.vars*dim(model.matrix(object))[2])
  cat("\n\n")
  cat("Residual Deviance:\t   ", format(signif(-2*sum(object$logL), digits)))
  cat("\n\n")
  cat("Dunn-Smyth Residuals:")
  cat("\n")
  cat(summary(qnorm(as.vector(object$residuals))))
  cat("\n\n-")
}
  
  
logLik.manyany <- function(object, ...)
{
  object$logL
}

residuals.manyany<- function(object, ...)
{
  tol=1.e-8
  params = object$params
  n.rows = length(params[[1]]$q)
  n.vars = length(params)
  if(length(object$family)==1)
    family = rep(object$family,n.vars)
  else
    family=object$family
  resids=matrix(NA,n.rows,n.vars)
  dimnames(resids)[[1]] = names(params[[1]]$yMat)
  dimnames(resids)[[2]] = names(params)
  for(i.var in 1:n.vars)
  {
    param.minus = params[[i.var]]
    param.minus$q = params[[i.var]]$q - 1.e-6
    if(family[[i.var]]$family=="negative.binomial")
      pfn = "pnbinom"
    if(family[[i.var]]$family=="poisson")
      pfn = "ppois"
    if(family[[i.var]]$family=="binomial")
      pfn = "pbinom"
    if(family[[i.var]]$family=="gaussian")
    {
      pfn = "pnorm"
      param.minus$q = params[[i.var]]$q
    } 
    if(family[[i.var]]$family=="Tweedie")
      pfn = "ptweedie"
    u = runif(n.rows)
    #to avoid any values identically 1:
    pMinus = pmin(do.call(pfn, param.minus), 1-tol)
    resids[,i.var] = u*do.call(pfn, params[[i.var]]) + (1-u)*pMinus
  }
  if(object$composition==TRUE) #reshape to original data size if required
    resids = matrix(resids, dim(object$fitted)[1], dim(object$fitted)[2])
  return(resids)
}


plot.manyany=function(x, ...)
{
  object = x
  Dunn.Smyth.Residuals=qnorm(residuals.manyany(object))
  Fitted.values=object$linear

  # add colours if not already there...
  if(hasArg("col")==F)
  {
    n.rows = dim(object$fitt)[1]
    n.vars = dim(object$fitt)[2]
    col=rep(1:n.vars,each=n.rows)
    plot(Dunn.Smyth.Residuals~Fitted.values, ann=F, col=col, ...)
  }
  else
    plot(Dunn.Smyth.Residuals~Fitted.values, ann=F, ...)
  
  #add xlab and ylab if not already provided...
  args=match.call()
  mx=match("xlab",names(args),0L)
  if(mx==0)
    xlab="Linear Predictor"
  else
    xlab=as.character(args[mx])

  my=match("ylab",names(args),0L)
  if(my==0)
    ylab="Dunn-Smyth Residuals"
  else
    ylab=as.character(args[my])
  
  title(xlab=xlab,ylab=ylab)

} 



