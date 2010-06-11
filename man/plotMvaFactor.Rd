\name{plotMvaFactor}
\alias{plotMvaFactor}
\title{Draw a Mvabund Object split into groups. }
\description{Draw the \code{mvabund} object \code{x} but split the data into
groups according to the grouping variable \code{y}. See \code{\link{default.plotMvaFactor}}
for a full description of the function.}
\usage{
plotMvaFactor(x, y, type="p", main="Abundance", n.vars= min(12,NCOL(x)),
	transformation="log", legend=TRUE, ...)
}
\arguments{ 
\item{x}{a \code{mvabund object}, a matrix with multivariate abundance data.}
\item{y}{a factor or a data.frame with factors, non-factor columns in a
data.frame are ignored.}
\item{type}{what type of plot should be drawn, allowed types are "p" for
scatterplot, "bx" for boxplot and "n" for no plot. Other types, as used
in \code{par} are NOT allowed.}
\item{main}{the title of the plot, see \code{plot}.}
\item{n.vars}{the number of variables to include in the plot.}
\item{transformation}{an optional transformation, "no" = untransformed,
"sqrt"=square root transformed, "log" (default)=log(Y/min+1) transformed,
"sqrt4" =4th root transformed.}
\item{legend}{logical, whether a legend should be added to the plot.}
\item{\dots}{arguments to be passed to or from other methods.}

}
\details{
For each variable in y that is a factor, a plot is drawn. When boxplots are drawn
the colors, that can be supplied by \code{col} are used to display different
factor levels. For scatterplots (type="p") it is also possible to use any plotting symbols
via the \code{pch} functional argument. \cr

If the colors and for scatterplots the plotting symbols are not supplied,
they will be automatically generated. If colors or the plotting symbols are supplied and 
if the number of factor levels is bigger than the the number of different values 
supplied by \code{pch} or \code{col}, then the colours and symbols will be repeated
for the factors. NOTE: This is also the default behaviour for the \code{plot} function.

Sometimes the legends might be only partially visible, especially when the width
of the graphics device is too small. To fix this, create a graphics device with
a larger width (see help("device") for available devices and their details)
and then repeat the \code{plotMvaFactor} command.
}
% \value{...}
\references{
Warton, D. I. (2008)
\emph{Raw data graphing: an informative but under-utilised tool 
for the analysis of multivariate abundances}, Austral Ecology.
}
% \note{...}
\author{
Ulrike Naumann, Stephen Wright and David Warton <David.Warton@unsw.edu.au>.
}
\seealso{
\code{\link{plot.mvabund}}.
}
\examples{
require(graphics)

#### Basic use ####
data(spider)
spiddat <- spider$abund
X <- spider$x

#Create a simplistic factor
X <- as.factor(X[,2]>0)
spiddat <- mvabund(spiddat)

plotMvaFactor(spiddat,X)		#complete default
plotMvaFactor(spiddat,X,type="bx")	#box plot via factor (TRUE or FALSE)

#### Advanced Use ####
data(solberg)
solbdat <- mvabund(solberg$abund)
treatment<- solberg$x

#Scatter plot version
plotMvaFactor(solbdat,treatment)

#Complex boxplot version
plotMvaFactor(x=solbdat, y=treatment, type="bx",
	main="BoxPlot of The 4 Highest Abundant Species",
	n.vars=4, xlab="Abundance [sqrt scale]", ylab="",
	transformation="sqrt", t.lab="o", shift=TRUE)

}
\keyword{hplot}

