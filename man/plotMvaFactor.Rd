\name{plotMvaFactor}
\alias{plotMvaFactor}
\title{Draw a Mvabund Object split into groups. }
\description{Draw the \code{mvabund} object \code{x} but split the data into
groups according to the grouping variable \code{y}.}
\usage{
plotMvaFactor(x, y, type="p", main="Abundances", xlab="Abundances", 
  ylab="Species", col= if(type=="bx") "white" else "black", 
  fg= "grey", pch=1, las=1, write.plot="show", filename="plot.mvabund", 
  n.vars= min(12,NCOL(x)), overall.main, var.subset=NA, subset=NA, 
  transformation="log", scale.lab="s", t.lab="o", mfrow=min(4,pExpl), 
  mfcol=NULL, shift=TRUE, border="black", keep.window=FALSE, ask=TRUE,
  legend=TRUE, legend.horiz=FALSE, legend.title=NULL, \dots)
}
\arguments{
\item{x}{a \code{mvabund object}, a matrix with multivariate abundance data.}
\item{y}{a factor or a data.frame with factors, non-factor columns in a
data.frame are ignored.}
\item{type}{what type of plot should be drawn, allowed types are "p" for
scatterplot, "bx" for boxplot and "n" for no plot. Other types, as used
in \code{par} are NOT allowed.}
\item{main}{the title of the plot, see \code{plot}.}
\item{xlab}{a title for the x axis, see \code{plot}, default is "Abundances".}
\item{ylab}{a title for the y axis, see \code{plot}, default is "Species".}
\item{col}{the color of the plotting symbols.
\code{col} might be used to make the different factor levels of the independent
variables visible. See \code{par} for a description of how colors are
specified and see the details for how to use colors to show the factor levels.}
\item{fg}{the color of the axis, default=grey.}
\item{pch}{the plotting symbol. See \code{points} for possible values and
their interpretation.
\code{pch} might be used to make the different factor levels of the independent
variables visible, see Details.}
\item{las}{the style of axis labels, see \code{par}.}
\item{write.plot}{
can be "show" (default): plot is shown, or "eps", "postscript", "pdf", "jpeg",
"bmp", "png" to save the plot in this format, plot is then not shown.}
\item{filename}{character, pathname where the file should be stored. 
Default location is R's working directory.}
\item{n.vars}{the number of variables to include in the plot.}
\item{overall.main}{a character to display as title for every window.}
\item{var.subset}{a numeric vector of indices indicating which variables of the
mvabund.object should be included on the plot.}
\item{subset}{a vector specifying a subset of cases. If any values in subset
are NA, the argument is discarded. Double variables are allowed.}
\item{transformation}{an optional transformation, "no" = untransformed,
"sqrt"=square root transformed, "log" (default)=log(Y/min+1) transformed,
"sqrt4" =4th root transformed.}
\item{scale.lab}{argument for limits and labelling of tickmarks:
"r" = R's default limits are used
"s" (default) = standard limits beginning at 0 are used.}
\item{t.lab}{argument for limits and labelling of tickmarks:
"o" (default) = "original" values of y are visible in the labels, e.g.
for transformation="sqrt", then the tickmark labels are based on x (and y)
not sqrt(x) (and sqrt(y))
"t" = transformed  values of y are visible in the labels.}
\item{mfrow}{layout of the plotting window, either a pair of c(rows, columns)
per window or an integer giving the prefered number of plots per window.}
\item{mfcol}{like \code{mfrow}, but the order of drawing is by row and not
by column, if mfcol is passed, mfrow is ignored.}
\item{shift}{logical, whether overlapping points in the plot should be shifted,
so that every point is visible.}
\item{border}{an optional value of colors for the outlines of the boxplots.}
\item{keep.window}{ whether, the specifications of mfrow, mfcol and windowsize
are kept.}
\item{ ask }{logical. If \code{TRUE} the user is asked for input,
before a new figure is drawn.}
\item{legend}{logical, whether a legend should be added to the plot.}
\item{legend.horiz}{logical, if applicable, whether the legend columns should
have horizontal orientation.}
\item{legend.title}{if applicable, the title of the legend.}
\item{\dots}{arguments to be passed to or from other methods.}
}
\details{
For each variable in y that is a factor, a plot is drawn. When boxplots are drawn
the colors, that can be supplied by \code{col} are used to display different
factor levels.
For scatterplots it is also possible to use the plotting symbols, specified by
\code{pch} for that. \cr

If the colors and for scatterplots the plotting symbols are not supplied,
they will be automatically generated. However, the plotting symbols will only
be automatically used in this way if there are up to seven different levels.

If colors or the plotting symbols are supplied, but the number of factor levels
is bigger than the the number of different values, they will be replicated.

Sometimes the legends might be only partially visible, especially when the width
of the graphics device is too small. To fix this, create a graphics device with
a larger width (see help("device") for on available devices and their details)
and then repeat the
\code{plotMvaFactor} command.
}
% \value{...}
\references{
Warton, D. I. ( )
\emph{Raw data graphing: an informative but under-utilised tool 
for the analysis of multivariate abundances}, , .
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

## Plot an Environment Factor vs Abundance plot
data(spider)
spiddat <- mvabund(spider$abund)
X <- spider$x

## Create a Environmental factor where TRUE=Sand, FALSE=No Sand)
X <- as.factor(X[,2]>0)
plotMvaFactor(x=spiddat, y=X)
}
\keyword{hplot}

