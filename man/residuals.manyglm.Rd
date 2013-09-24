\name{residuals.manyglm}
\alias{residuals.manyglm}
\title{Residuals for MANYGLM Fits}
\description{  
Obtains Dunn-Smyth residuals from a fitted manyglm object.
}
\usage{
% residuals.manyglm
\method{residuals}{manyglm}(object, \dots)
}
\arguments{
  \item{object}{a fitted object of class inheriting from \code{"manyglm"}.}
  \item{\dots}{further arguments passed to or from other methods.}
}
\details{
residuals.manyglm computes Randomised Quantile or ``Dunn-Smyth" residuals (Dunn
& Smyth 1996) for a manyglm object. If the fitted model is correct then Dunn-Smyth residuals
are standard normal in distribution.

Note that for discrete data, Dunn-Smyth residuals involve random number generation, and 
will not return identical results on replicate runs. Hence it is worth calling this function
multiple times to get a sense for whether your interpretation of results holds up under replication.
}
\value{
  A matrix of Dunn-Smyth residuals.  
}
\references{

Dunn, P.K., & Smyth, G.K. (1996). Randomized quantile residuals. Journal of Computational and
Graphical Statistics 5, 236-244.

}
\author{David Warton <David.Warton@unsw.edu.au>.}
\seealso{
\code{\link{manyglm}}.
}
\examples{
data(spider)
spiddat <- mvabund(spider$abund)
X <- spider$x

## obtain residuals for Poisson regression of the spider data, and doing a qqplot:
glmP.spid  <- manyglm(spiddat~X, family="poisson")
resP       <- residuals(glmP.spid)
qqnorm(resP)
qqline(resP,col="red")
#clear departure from normality.

## try again using negative binomial regression:
glmNB.spid <- manyglm(spiddat~X, family="negative.binomial")
resNB      <- residuals(glmNB.spid)
qqnorm(resNB)
qqline(resNB,col="red")
#that looks a lot more promising.

#note that you could construct a similar plot directly from the manyglm object using
plot(glmNB.spid, which=2)

}
\keyword{models}
\keyword{multivariate}
\keyword{regression}

