\name{mvabund-package}

\alias{mvabund-package}

\title{Statistical methods for analysing multivariate abundance data }

\description{
This package provides tools for a model-based approach to the analysis of multivariate abundance data in ecology (Warton (2011)).

This includes graphical methods for exploring the properties of data and the community-environment association, flexible regression methods for estimating and making robust inferences about the community-environment association, and diagnostic plots to check the appropriateness of a fitted model.

Future versions of this function will extend it to regression analysis using additive functions, and model-based cluster analysis.
}

\details{
The key functions available in this package are the following.

\bold{For graphical display of the data:}
    \describe{
      \item{\code{\link{plot.mvabund}}}{ draw a range of plots for Multivariate Abundance Data }
      \item{\code{\link{boxplot.mvabund}}}{ draw a range of plots of Model Formulae for Multivariate Abundance Data }
      \item{\code{\link{meanvar.plot}}}{ draw mean-variance plots for Multivariate Abundance Data }
      \item{\code{\link{plotMvaFactor}}}{ draw a range of plots for Multivariate Abundance Data where the data is split into groups according to a grouping variable}
    }

\bold{For estimating and displaying Linear Models:}
    \describe{
      \item{\code{\link{manylm}}}{Fitting Linear Models for Multivariate Abundance Data}
      \item{\code{\link{summary.manylm}}}{summarizie Multivariate Linear Model Fits for Abundance Data}
      \item{\code{\link{anova.manylm}}}{obtain ANOVA for Multivariate Linear Model Fits for Abundance Data}
      \item{\code{\link{plot.manylm}}}{plot diagnostics for a manylm Object}
   }

\bold{For estimating and displaying Generalized Linear Models:}
    \describe{
      \item{\code{\link{manyglm}}}{fit Generalized Linear Models for Multivariate Abundance Data}
      \item{\code{\link{summary.manyglm}}}{summarize Multivariate Generalized Linear Model Fits for Abundance Data}
      \item{\code{\link{anova.manyglm}}}{obtain Analysis of Deviance for Multivariate Generalized Linear Model Fits for Abundance Data}
      \item{\code{\link{plot.manyglm}}}{plot diagnostics for a manyglm Object}
    }

\bold{For providing a data structure:}
    \describe{
      \item{\code{\link{mvabund}}}{create a mvabund object}
      \item{\code{\link{mvformula}}}{create Model Formulae for Multivariate Abundance Data}
    }
  

\bold{Example datasets:}
    \describe{
      \item{\code{\link{Tasmania}}}{ meiobenthic community data from Tasmania. Used to demonstrate test for interaction. }
      \item{\code{\link{solberg}}}{ solberg species data with a 3-level treatment factor. }
      \item{\code{\link{spider}}}{ hunting spiders data from different sites.}
      \item{\code{\link{tikus}}}{ solberg nematode data from Tikus island.}
    }

For more details, see the documentation for any of the individual functions listed above.

}
\references{

Lawless, J. F. (1987) \emph{Negative binomial and mixed Poisson regression}, Canadian Journal of Statistics 15, 209-225.

Liang, K.-Y. and Zeger, S. L (1986) \emph{Longitudinal Data Analysis Using Generalized Linear Models}, Biometrika 73, No. 1, 13-22.

Warton D.I. (2008a). Raw data graphing: an informative but under-utilized tool for the analysis of multivariate abundances. \emph{Austral Ecology} 33(3), 290-300.

Warton D.I. (2008b). Penalized normal likelihood and ridge regularization of correlation and covariance matrices. \emph{Journal of the American Statistical Association} 103, 340-349.

Warton D.I. & Hudson H.M. (2004). A MANOVA statistic is just as powerful as distance-based statistics, for multivariate abundances. \emph{Ecology} 85(3), 858-874.

Warton, D. I. and Guttorp, P. (in press April 2010). Compositional analysis of overdispersed counts using generalized estimating equations. \emph{Environmental and Ecological Statistics}.

Warton D.I. (2011). Regularized sandwich estimators for analysis of high dimensional data usi
ng generalized estimating equations. \emph{Biometrics}, 67(1), 116-123.

}

\author{David Warton \email{David.Warton@unsw.edu.au}, Yi Wang and Ulrike Naumann.}

\seealso{ 
\code{\link{mvabund}}, \code{\link{mvformula}},
\code{\link{plot.mvabund}}, \code{\link{boxplot.mvabund}}, \code{\link{meanvar.plot}}, \code{\link{plotMvaFactor}},
\code{\link{manylm}}, \code{\link{summary.manylm}}, \code{\link{anova.manylm}}, \code{\link{plot.manylm}}, 
}

\examples{
require(graphics)

## Load the spider dataset:
data(spider)

## Create the mvabund object spiddat:
spiddat <- mvabund(spider$abund)
X <- spider$x

## Draw a plot of the spider data:
plot(spiddat, col="gray1", n.vars=8, transformation="sqrt", xlab=c("Hunting Spider"), ylab="Spider Species", scale.lab="s",t.lab="t", shift=TRUE, fg= "lightblue", col.main="red", main="Spiders") 


## A mean-variance plot, data organised by year, 
## for 1981 and 1983 only, as in Figure 7a of Warton (2008):
data(tikus)
tikusdat <- mvabund(tikus$abund)
year <- tikus$x[,1]
is81or83 <- year==81 | year==83
meanvar.plot(tikusdat~year,legend=TRUE, subset=is81or83, col=c(1,10)) 	

## Create a formula for multivariate abundance data:
foo <- mvformula( spiddat~X )

## Create a List of Univariate Formulas:
fooUni <- formulaUnimva(spiddat~X)
fooUniInt <- formulaUnimva(spiddat~X, intercept=TRUE)

## Find the three variables that best explain the response:
best.r.sq( foo, n.xvars= 3)

## Fit a multivariate linear model:
foo <- mvformula( spiddat~X )
lm.spider <- manylm(foo)

## Plot Diagnostics for a multivariate linear model:
plot(lm.spider,which=1:2,col.main="red",cex=3,overlay=FALSE)

## Obtain a summary of test statistics using residual resampling:
summary(lm.spider, nBoot=500)

## Calculate a ANOVA Table:
anova(lm.spider, nBoot=500)

}
    
\docType{ package }
\keyword{ documentation }
