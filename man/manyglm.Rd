\name{manyglm}
\alias{manyglm}
\title{Fitting Generalized Linear Models for Multivariate Abundance Data}
\description{
  \code{manyglm} is used to fit generalized linear models to high-dimensional data, such as multivariate abundance data in ecology. This is the base model-fitting function - see \code{plot.manyglm} for assumption checking, and \code{anova.manyglm} or \code{summary.manyglm} for significance testing.
}
\usage{
manyglm(formula, family="negative.binomial", K=1, data=NULL, subset=NULL, 
    na.action=options("na.action"), theta.method = "PHI", model = FALSE,
    x = TRUE, y = TRUE, qr = TRUE, cor.type= "I", shrink.param=NULL, 
    tol=sqrt(.Machine$double.eps), maxiter=25, maxiter2=10, show.coef=FALSE, show.fitted=FALSE, show.residuals=FALSE, show.warning=FALSE,\dots )
}
\arguments{
  \item{formula}{an object of class \code{"\link{formula}"} (or one that
    can be coerced to that class): a symbolic description of the
    model to be fitted.  The details of model specification are given
    under Details.}
  \item{family}{a description of the error distribution function to be used 
    in the model. At the moment, this must be a character string naming one 
    of the following members of the exponential family: 'gaussian', 'poisson',
    'binomial', 'negative.binomial'. Further improvement will include options
    of specifying the link function and calling to a family function defined 
    in \code{\link{family}}.
   }
\item{K}{number of trials in binomial regression. By default, K=1 for presence-absence data using logistic regression.}
   \item{data}{an optional data frame, list or environment (or object
    coercible by \code{\link{as.data.frame}} to a data frame) containing
    the variables in the model. If not found in \code{data}, the variables 
    are taken from \code{environment(formula)}, typically the environment 
    from which \code{glm} is called.}
  \item{subset}{an optional vector specifying a subset of observations
    to be used in the fitting process.}
  \item{na.action}{a function which indicates what should happen
    when the data contain \code{NA}s.  The default is set by
    the \code{na.action} setting of \code{\link{options}}, and is
    \code{\link{na.fail}} if that is unset.  The \sQuote{factory-fresh}
    default is \code{\link{na.omit}}.  Another possible value is
    \code{NULL}, no action.  Value \code{\link{na.exclude}} can be useful.}
  \item{theta.method}{the method used for the estimation of the overdisperson 
  parameter theta, such that the mean-variance relationship is V=m+m^2/theta for 
  the negative binomial family. Here offers two options \cr
  "ML" = Maximum likelihood estimation (default), as in Lawless(1987), \cr
  "Chi2" = Moment estimation using chi-square dampening on the log scale, as
  in Hilbe(2008).
  %"MLFisher" (= ML with Fisher scoring, only available for the negative.binomial
  %family, the most efficient of the ML methods) and
  % "MLLawless" (ML using iterative estimation as in Lawless (1987),
  % only available for negative.binomial family).
  }
  \item{model, x, y, qr}{logicals. If \code{TRUE} the corresponding 
  components of the fit (the model frame, the model matrix, the model
  matrix, the response, the QR decomposition of the model matrix) are
  returned. 
  }
%  \item{test}{the test statistic to be used in subsequent functions for
%  inference. It is merely stored in \code{manyglm} as an attribute. 
%  Possible options are\cr
%  "LR" = Log-likelihood ratio test (default),\cr
%  "F" = Lawley-Hotelling trace statistic for 'gaussian' only,\cr
%  "wald" = Wald statistic (not available for 'gaussian' yet), \cr
%  "score" = Score statistic (not available for 'gaussian' yet).
%  }
  \item{cor.type}{the structure imposed on the estimated correlation 
  matrix under the fitted model. Can be "I"(default), "shrink", or "R". 
  See Details. This parameter is merely stored in \code{manyglm}, and 
  will be used as the default value for \code{cor.type} in subsequent 
  functions for inference.
  }
  \item{shrink.param}{
  shrinkage parameter to be used if \code{cor.type="shrink"}. If a numerical 
  value is not supplied, it will be estimated from the data by cross 
  validation-penalised normal likelihood as in Warton (2008). The parameter 
  value is stored as an attribute of the \code{manyglm} object, and will be 
  used in subsequent functions for inference.
  }
  \item{tol}{the tolerance used in estimations.}
  \item{maxiter}{maximum allowed iterations in the weighted least square estimation of beta. The default value is 100.}
  \item{maxiter2}{maximum allowed iterations in the internal ML estimations of negative binomial regression. The default value is 5.}
  \item{show.coef, show.fitted, show.residuals, show.warning}{logical. Whether to show model coefficients, fitted values, standardized pearson residuals, or operation warnings.}
  \item{\dots}{further arguments passed to or from other methods.}
}
\details{

\code{manyglm} is used to calculate the parameter estimates of generalised linear models fitted to each of many variables simultaneously as in Warton et. al. (2012) and Wang et.al.(2012). Models for \code{manyglm} are specified symbolically. For details on this compare the details section of \code{\link{lm}} and \code{\link{formula}}. A formula has an implied intercept term.  To remove this use either \code{y ~ x - 1} or \code{y ~ 0 + x}.  See \code{\link{formula}} for more details of allowed formulae.

Generalised linear models are designed for non-normal data where a specific alternative distribution offers a reasonable model for data, as specified using the argument \code{family} (Warton (2011)). Count data can be analysed using \code{family="poisson"} if the data is not "overdispersed" (that is, if the variance is not larger than the mean), although for multivariate abundance data it has been shown that the negative binomial distribution (\code{family="negative.binomial"}) is a good choice of model for counts (Warton (2005)). For binary (presence/absence) data, \code{family="binomial"} should be used.

The negative binomial function requires knowledge of the overdispersion
parameter (\code{theta}) in order to fit a generalised linear model. However, it
is usual not to have knowledge of the value of \code{theta}. So we estimate 
%so if it is not specified it will be estimated from the data. The default estimates 
a separate \code{theta} for each variable from the data.
The method of estimating \code{theta} is controlled by \code{theta.method} for negative binomial distributions. 

\code{family} not only chooses the type of distribution to fit to data,
but also controls the link function - the function of the mean that is fitted
with a linear predictor. The default is the canonical link function, except
for negative binomial, which uses the log-link. 
%The variance-stabilising parameterisation (\code{link="varstab"}) is a special type of link function that was introduced to the \code{manyglm} function because it allows construction of Wald statistics with good properties (Warton 2008), in situations where its application is reasonable, such as multi-sample testing.

%If \code{family=quasipoisson} is used (that is, the mean-variance relationship
%is V=phi*m), \code{phi} can be estimated from data as for
%\code{family=negative.binomial}, as well as specified in advance. For
%\code{family=quasibinomial} however \code{phi} cannot be specified in advance.

%Note that when a family \code{binomial} or \code{quasibinomial} is specified in \code{manyglm}, the workhorse-functions used for the family are actually \code{binomialMany} and \code{quasibinomialMany} respectively.

%The response can be in the one of two formats: \cr
%Either the response is given in a matrix of bernoulli variables, that is each
%element of the matrix either has the value one or zero. \cr
%The response may also be given as matrix that combines a matrix with numbers of successes and a matrix with number of failures. This matrix can be produced by using \code{\link{binstruc}}. The failures and successes must be identifiable by the column names, which start either with the string \dQuote{succ} or \dQuote{fail}. The order of variables must be the same in the matrix of successes and the matrix of failures and the number of trials (i.e. the sum of failures and successes) must be the same for corresponding rows of each variable.
% (This may change in further version of the \code{mvabund} package).

% A terms specification of the form \code{first + second} indicates all the terms in \code{first} together with all the terms in \code{second} with duplicates removed. The terms in the formula will be re-ordered so that main effects come first, followed by the interactions, all second-order, all third-order and so on: to avoid this pass a \code{terms} object as the formula.

%Then \code{manyglm.fit} is the workhorse function. For other families the default is to using "glm.fit". Then the model is rewritten as a univariate model and the function 'glm.fit' will be used to do the actual computation. Depending on the dimension of the data this can be faster than 'manyglm.fit', however cannot be used for "poisson", "quasipoisson" and "Negative Binomial". Both methods use iteratively reweighted least squares (IWLS).

%If more than one of \code{etastart}, \code{start} and \code{mustart}
%is specified, the first in the list will be used.  It is often
%advisable to supply starting values for a \code{\link{quasi}} family,
%and also for families with unusual links such as \code{gaussian("log")}.

% when 'manyglm.fit' is used none of \code{etastart}, \code{start} and \code{mustart} is actually used except in following testing functions (summary.manyglm, ...)

%All of \code{weights}, \code{subset}, \code{offset}, \code{etastart} and
%\code{mustart} are evaluated in the same way as variables in \code{formula},
%that is first in \code{data} and then in the environment of \code{formula}.

\code{cor.type} is the structure imposed on the estimated correlation
matrix under the fitted model. Possible values are: \cr
"I"(default) = independence is assumed (correlation matrix is the identity) \cr
"shrink" = sample correlation matrix is shrunk towards I to improve its stability. \cr
"R" = unstructured correlation matrix is used.  (Only available when N>p.)

If \code{cor.type=="shrink"}, a numerical value will be assigned to \code{shrink.param} either through the argument or by internal estimation. The working horse function for the internal estimation is \code{\link{ridgeParamEst}}, which is based on cross-validation (Warton 2008). The validation groups are chosen by random assignment, so some slight variation in the estimated values may be observed in repeat analyses. See \code{\link{ridgeParamEst}} for more details. The shrinkage parameter can be any value between 0 and 1 (0="I" and 1="R", values closer towards 0 indicate more shrinkage towards "I").
}

\value{
  \code{manyglm} returns an object inheriting from \code{"manyglm"},
  \code{"manylm"} and "mglm".

  The function \code{summary} (i.e. \code{\link[mvabund:anova.manyglm]{summary.manyglm}}) can be used to obtain or print a summary of the results and the function
  \code{anova} (i.e. \code{\link[mvabund:anova.manyglm]{anova.manyglm}}) to produce an 
  analysis of variance table.

  The generic accessor functions \code{\link{coefficients}}, 
  \code{fitted.values} and \code{residuals} can be used to
  extract various useful features of the value returned by \code{manyglm}.

% \code{weights} extracts a vector of weights, one for each case in the fit (after subsetting and \code{na.action}).

  An object of class \code{"manyglm"} is a list containing at least the
  following components:

  \item{coefficients}{a named matrix of coefficients.}
  \item{fitted.values}{the matrix of fitted mean values, obtained by transforming the linear predictors by the inverse of the link function.}
  \item{linear.predictors}{the linear fit on link scale.}
  \item{residuals}{the matrix of \emph{working} residuals, that is the Pearson residuals standardized by the leverage adjustment h obtained from the diagonal elements of the hat matrix H. }
  \item{sqrt.1_Hii}{the matrix of scale terms used to standardize the Pearson reidusals.}
  \item{variance.estimator}{the variance estimator for the corresponding family function.}
  \item{sqrt.weight}{the matrix of square root of \emph{working} weights,  estimated for the corresponding family function.}
  \item{theta}{the estimated nuisance parameters accounting for overdispersion}
  \item{two.loglike}{two times the log likelihood.}
  \item{deviance}{up to a constant, minus twice the maximized log-likelihood.  Where sensible, the constant is chosen so that a saturated model has deviance zero.}
  \item{aic}{Akaike's \emph{An Information Criterion}, minus twice the
    maximized log-likelihood plus twice the number of coefficients.} 
  \item{AICsum}{the sum of the AIC's over all variables.}
  \item{iter}{the number of iterations of IWLS used.}
  \item{tol}{the tolerance used in estimations.}

  \item{family}{the \code{\link{family}} argument supplied.}
  \item{theta.method}{the \code{theta.method} argument supplied.}
  \item{cor.type}{the \code{cor.type} argument supplied.}
  \item{shrink.param}{the shrink parameter to be used in subsequent inference.}
  \item{call}{the matched call.}
  \item{terms}{the \code{\link{terms}} object used.}
  \item{formula}{the formula supplied.}
  \item{rank}{the numeric rank of the fitted linear model.}
  \item{xlevels}{(where relevant) a record of the levels of the factors
    used in fitting.}
  \item{df.residual}{the residual degrees of freedom.}

  \item{model}{if the argument \code{model} is \code{TRUE}, this is the \code{model.frame}.}
  \item{y}{if the argument \code{y} is \code{TRUE}, this is the response variables used.}
  \item{x}{if the argument \code{x} is \code{TRUE}, this is the design matrix used.}
  \item{qr}{if the argument \code{qr} is \code{TRUE}, this is the QR decomposition of the design matrix.}  

% In addition, non-empty fits will have components \code{qr} and \code{R}, both
%  being lists with components corresponding to the dependent variables.
%  Furthermore, non-empty fits will have \code{effects} relating to the final
%  weighted linear fit.

  Objects of class \code{"glm"} are normally of class \code{c("glm",
    "lm")}, that is inherit from class \code{"lm"}, and well-designed
  methods for class \code{"lm"} will be applied to the weighted linear
  model at the final iteration of IWLS.  However, care is needed, as
  extractor functions for class \code{"glm"} such as
  \code{\link{residuals}} and \code{weights} do \bold{not} just pick out
  the component of the fit with the same name.

%  If a \code{manyglm} model with a binomial family was specified by giving a
%  two-column response, the weights returned by \code{prior.weights} are
%  the total numbers of cases (factored by the supplied case weights) and
%  the component \code{y} of the result is the proportion of successes.
}

\references{


Lawless, J. F. (1987)
\emph{Negative binomial and mixed Poisson regression},
Canadian Journal of Statistics 15, 209-225.

Hilbe, J. M. (2008)
\emph{Negative Binomial Regression},
Cambridge University Press, Cambridge.

Warton D.I. (2005)
\emph{Many zeros does not mean zero inflation: comparing the
goodness of-fit of parametric models to multivariate abundance data},
Environmetrics 16(3), 275-289.

Warton D.I. (2008). Penalized normal likelihood and ridge regularization
of correlation and covariance matrices. \emph{Journal of the American
Statistical Association} 103, 340-349.

Warton D.I. (2011). Regularized sandwich estimators for analysis of high dimensional data using generalized estimating equations. \emph{Biometrics}, 67(1), 116-123.

Warton D. I., Wright S., and Wang, Y. (2012). Distance-based multivariate analyses confound location and dispersion effects. \emph{Methods in Ecology and Evolution}, 3(1), 89-101.

Wang Y., Neuman U., Wright S. and Warton D. I. (2012). mvabund: an R package for model-based analysis of multivariate abundance data. \emph{Methods in Ecology and Evolution}. online 21 Feb 2012.

}

\author{
Yi Wang, Ulrike Naumann and David Warton <David.Warton@unsw.edu.au>.
}
\seealso{
%\code{\link{negative.binomial}}, \code{\link{quasipoisson}},
%\code{\link{binomialMany}}, \code{\link{quasibinomialMany}},
%\code{\link{binstruc}}, 
\code{\link[mvabund:anova.manyglm]{anova.manyglm}}, \code{\link[mvabund:summary.manyglm]{summary.manyglm}},
%\code{\link{fitted.values}}, \code{\link{residuals}}.
}
\examples{
data(spider)
spiddat <- mvabund(spider$abund)
X <- spider$x

#To fit a log-linear model assuming counts are poisson:
glm.spid <- manyglm(spiddat~X, family="poisson")
glm.spid 

summary(glm.spid, resamp="residual")

#To fit a binomial regression model to presence/absence data:
pres.abs <- spiddat
pres.abs[pres.abs>0] = 1
glm.spid.bin <- manyglm(pres.abs~X, family="binomial")
glm.spid.bin

summary(glm.spid.bin, nBoot=500)


}
\keyword{models}
\keyword{regression}
\keyword{multivariate}

