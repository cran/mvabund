\name{anova.manyglm}
\alias{anova.manyglm}
\alias{print.anova.manyglm}
\title{Analysis of Deviance for Multivariate Generalized Linear Model
Fits for Abundance Data}
\description{Compute an analysis of deviance table for one or more multivariate
generalized linear model fits.}
\usage{
\method{anova}{manyglm}(object, \dots, resamp="residual", test="LR", 
    p.uni="none", nBoot=1000, cor.type=object$cor.type, ld.perm=FALSE, filename=NULL) 
\method{print}{anova.manyglm}(x, \dots) 
}
\arguments{
\item{object}{objects of class \code{manyglm}, typically the result of a call to \code{\link{manyglm}}.
}
\item{\dots}{for the \code{anova.manyglm} method, these are optional further objects of class \code{manyglm}, which are usually a result of a call to \code{\link{manyglm}}
for the \code{print.anova.manyglm} method these are optional further arguments passed to or from other methods. See \code{\link{print.summary.glm}} for more details.}        
\item{resamp}{the method of resampling used. Can be one of "case", "residual" (default), "perm.resid", "score", or "montecarlo". See Details.}
\item{test}{the test to be used. If \code{test="I"}, this can be one of "wald" for a Wald-Test or "score" for a Score-Test or "LR" for a Likelihood-Ratio-Test, otherwise only "wald" and "score" is allowed. The default value is "LR".
}
\item{p.uni}{whether to calculate univariate test statistics and their P-values, and if so, what type. This can be one of the following options. \cr
"none" = No univariate P-values (default) \cr
"unadjusted" = A test statistic and (ordinary unadjusted) P-value are reported
for each response variable. \cr
"adjusted" = Univariate P-values are adjusted for multiple testing, using a
step-down resampling procedure.
}
\item{nBoot}{the number of Bootstrap iterations, default is \code{nBoot=999}.}
\item{cor.type}{structure imposed on the estimated correlation matrix under the fitted model. Can be "I", "shrink", or "R". By default, the value is inherited from the first \code{\link{manyglm}} object. See Details.
}
\item{ld.perm}{logical. Whether to load bootstrap id's from an external file. Default is \code{FALSE}, which will generate bootstrap id's on the fly.}
\item{filename}{if \code{ld.perm} is \code{TRUE}, this specifies the file path to the stored bootstrap id's. Default is \code{NULL}.}
\item{x}{an object of class "summary.manyglm", usually, a result of a call to
\code{\link[mvabund:summary.manyglm]{summary.manyglm}}.}
}

\details{
The \code{anova.manyglm} function returns a table summarising the statistical significance of a fitted manyglm model (Warton (2011)), or of the differences between several nested models. If one model is specified, sequential test statistics (and P values) are returned for that fit. If more than one object is specified, the table contains test statistics (and P values) comparing their fits, provided that the models are fitted to the same dataset.

The test statistics are determined by the argument \code{test}, and the P-values are calculated by resampling rows of the data using a method determined by the argument \code{resamp}. The five possible resampling methods (case, residual, score resampling, residual permutation and parametric bootstrap) are described in more detail in Davison and Hinkley (1997, chapter 6) and Wu (1986), and involve resampling under the null hypothesis (except for case resampling). These methods ensure approximately valid inference even when the mean-variance relationship or the correlation between variables has been misspecified. Standardized pearson residuals (see \code{\link{manyglm}} are used in residual and score resampling, and where necessary, resampled response values are truncated so that they fall in the required range (e.g. counts cannot be negative). If \code{resamp="none"}, p-values cannot be calculated, however the test statistics are returned. For abundance data with a very high dimension and lots of zero values it is not recommended to use \code{residual} resampling with an \code{LR} test.

If you do not have a specific hypothesis of primary interest that you want to test, and are instead interested in which model terms are statistically significant, then the \code{summary.manyglm} function is more appropriate. Whereas \code{summary.manyglm} tests the significance of each explanatory variable, \code{anova.manyglm}, given one \code{manyglm} object tests each  term of the formula, e.g. if the formula is 'y~a+b' then a and b, that can be vectors or matrices, are tested for significance.

For information on the different types of data that can be modelled using manyglm, see \code{\link{manyglm}}. To check model assumptions, use \code{\link{plot.manyglm}}.

Multivariate test statistics are constructed using one of three methods: a log-likelihood ratio statistic \code{test="LR"} or a Wald statistic \code{test="wald"} or a Score statistic \code{test="score"}. "LR" has good properties, but is only available when \code{cor.type="I"}. 
%Also recommended, but only for multi-sample testing, is a Wald statistic using the variance-stabilising link function (e.g. \code{family=negative.binomial(varstab)}), which has been shown to have good power properties (Warton 2008b). However use of this link function is only reasonable for multi-sample tests, a special case where changing link function can be considered as a reparameterisation of the model (rather than a change of model). 

The default Wald test statistic makes use of a generalised estimating equations (GEE) approach, estimating the covariance matrix of parameter estimates using a sandwich-type estimator that assumes the mean-variance relationship in the data is correctly specified and that there is an unknown but constant correlation across all observations. Such assumptions allow the test statistic to account for correlation between variables but to do so in a more efficient way than traditional GEE sandwich estimators (Warton 2008a). The common correlation matrix is estimated from standardized Pearson residuals, and the method specified by \code{cor.type} is used to adjust for high dimensionality. 

The Wald statistic has problems for count data and presence-absence 
data when there are zero parameters, so is not recommended for multi-sample 
tests, where such situations are common. 

% By default, \code{anova.manyglm} tests hypotheses concerning mean (absolute) abundance, but it is also possible to test hypotheses concerning community composition using the argument \code{compositional=TRUE}, when analysing count data, as in Warton and Guttorp (2007). This constructs Wald statistics to test hypotheses concerning relative abundance not absolute abundance, and the parameters that are returned in the fitted model are for the relative change in abundance of one variable as compared to another (as compared to the last abundance variable in the mvabund array).

The \code{anova.manyglm} function is designed specifically for high-dimensional data (that, is when the number of variables p is not small compared to the number of observations N). In such instances a correlation matrix is computationally intensive to estimate and is numerically unstable, so by default the test statistic is calculated assuming independence of variables (\code{cor.type="I"}). Note however that the resampling scheme used ensures that the P-values are approximately correct even when the independence assumption is not satisfied. However if it is computationally feasible for your dataset, it is recommended that you use \code{cor.type="shrink"} to account for correlation between variables, or \code{cor.type="R"} when p is small. The \code{cor.type="R"} option uses the unstructured correlation matrix (only possible when N>p), such that the standard  classical multivariate test statistics are obtained. Note however that such statistics are typically numerically unstable and have low power when p is not small compared to N. 

The \code{cor.type="shrink"} option applies ridge regularisation (Warton 2008b), shrinking the sample correlation matrix towards the identity, which improves its stability when p is not small compared to N. This provides a compromise between \code{"R"} and \code{"I"}, allowing us to account for correlation between variables, while using a numerically stable test statistic that has good properties. 

The shrinkage parameter is an attribute of a \code{\link{manyglm}} object. For a Wald test, the sample correlation matrix of the alternative model is used to calculate the test statistics. So \code{shrink.param} of the alternative model is used. For a score test, the sample correlation matrix of the null model is used to calculate the test statistics. So \code{shrink.param} of the null model is used instead. If \code{cor.type=="shrink"} and \code{shrink.param} is NULL, then the shrinkage parameter will be estimated by cross-validation using the multivariate normal likelihood function (see \code{\link{ridgeParamEst}} and (Warton 2008b)) for the corresponding model in the anova test.

Rather than stopping after testing for multivariate effects, it is often of interest to find out which response variables express significant effects. Univariate statistics are required to answer this question, and these are reported if requested. Setting \code{p.uni="unadjusted"} returns resampling-based univariate P-values for all effects as well as the multivariate P-values, whereas \code{p.uni="adjusted"} returns adjusted P-values (that have been adjusted for multiple testing), calculated using a step-down resampling algorithm as in Westfall & Young (1993, Algorithm 2.8). This method provides strong control of family-wise error rates, and makes use of resampling (using the method controlled by \code{resamp}) to ensure inferences take into account correlation between variables.
}

\value{
\item{family}{the \code{family} component from \code{object}.}
\item{p.uni}{the \code{p.uni} argument supplied.}
\item{test}{the \code{test} argument supplied.}
\item{cor.type}{the \code{cor.type} argument supplied.}
\item{resamp}{the \code{resamp} argument supplied.}
\item{nBoot}{the \code{nBoot} argument supplied.}
\item{shrink.parameter}{a list of shrink parameters from all \code{manyglm} objects in the anova test.}
\item{n.bootsdone}{the number of bootstrapping iterations that were done, i.e. had no error.}
\item{table}{the table with Residual Degrees of Freedom, Degrees of Freedom, the Test Statistics and the P values.}

If \code{p.uni="adjusted"} or \code{"unadjusted"} the output list also contains\cr
\item{uni.test}{a table showing the test statistics of the univariate tests.}
\item{uni.p}{a table showing the p-values of the univariate tests.}
}

\section{Warning}{
The comparison between two or more models by \code{anova.manyglm} will only be valid if they are fitted to the same dataset. This may be a problem if there are missing values and R's default of \code{na.action = na.omit} is used. 
}

\references{
Warton D.I. (2011). Regularized sandwich estimators for analysis of high dimensional data using generalized estimating equations. \emph{Biometrics}, 67(1), 116-123.

Warton D.I. (2008a). Which Wald statistic?  Choosing a parameterisation of the Wald statistic to maximise power in \emph{k}-sample generalised estimating equations. \emph{Journal of Statistical Planning and Inference}, 138, 3269-3282.

Warton D.I. (2008b). Penalized normal likelihood and ridge regularization of correlation and covariance matrices. \emph{Journal of the American Statistical Association} 103, 340-349.

Warton D.I. and Guttorp P. (in press April 2010) Compositional analysis of overdispersed 
counts using generalized estimating equations. \emph{Environmental and Ecological Statistics} 
Davison, A. C. and Hinkley, D. V. (1997) \emph{Bootstrap Methods and their Application.} Cambridge University Press, Cambridge.

Westfall, P. H. and Young, S. S. (1993) \emph{Resampling-based multiple 
testing.} John Wiley & Sons, New York.

Wu, C. F. J. (1986) Jackknife, Bootstrap and Other Resampling Methods in
Regression Analysis. \emph{The Annals of Statistics} 14:4, 1261-1295.
}
\author{
Yi Wang, Ulrike Naumann and David Warton <David.Warton@unsw.edu.au>.
}
\seealso{
\code{\link{manyglm}}, \code{\link[mvabund:summary.manyglm]{summary.manyglm}}.
}
\examples{
## Load the Tasmania data set
data(Tasmania)

## Visualise the effect of treatment on copepod abundance
tasm.cop <- mvabund(Tasmania$copepods)
treatment <- Tasmania$treatment
block <- Tasmania$block
plot(tasm.cop ~ treatment, col=as.numeric(block))

## Fitting predictive models using a negative binomial model for counts:
tasm.cop.nb <- manyglm(tasm.cop ~ block*treatment, family="negative.binomial")

## Testing hypotheses about the treatment effect and treatment-by-block interactions, 
## using the default settings of residual resampling and deviance as test statistics: 
anova(tasm.cop.nb)

## Checking model assumptions
plot(tasm.cop.nb)


## Accounting for correlation between variables (with shrinkage):
tasm.cop.nb2 <- manyglm(tasm.cop ~ block*treatment, family="negative.binomial", cor.type="shrink")
## Repeat ANOVA test with P-values adjusted for multiple testing, using montecarlo parametric bootstrap and wald statistics:
anova(tasm.cop.nb2, p.uni="adjusted", resamp="montecarlo", test="wald")

}
\keyword{models}
\keyword{regression}
\keyword{multivariate}
