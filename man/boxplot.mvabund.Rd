\name{boxplot.mvabund}
\alias{boxplot.mvabund}
\alias{boxplot.mvformula}
\title{Boxplots for multivariate abundance Data}
\description{
Draw Boxplots of \code{mvabund} or \code{mvformula} Objects
}
\usage{
\method{boxplot}{mvabund}(
    x, y=NULL, ..., range=1.5, width=NULL, notch=FALSE, outline=TRUE, 
    names, border=par("fg"), col=NULL, pars=list(boxwex=0.8, 
    staplewex=0.5, outwex=0.5), at=NULL, xlab, ylab, main, pch=1, 
    fg="grey",  las=1, write.plot="show", filename="plot.mvabund", 
    n.vars=min(12,NCOL(x)), overall.main, var.subset=NA, 
    transformation="log", scale.lab="s", t.lab="o", mfrow=1, mfcol=NULL,
    ask=if(write.plot=="show"&(is.list(x)|!is.null(dotsnotnamed)))
    TRUE else FALSE )

\method{boxplot}{mvformula}(
    x, \dots, data=NULL, subset, main="", xlab=NULL, ylab="Abundances", 
    col="white", fg="grey", pch=20, las=NULL, n.vars=12, overall.main="", 
    var.subset=NA, write.plot="show", filename="plot.mvabund", 
    scale.lab=NULL, t.lab="o", mfrow=c(min(4,n.vars),3), mfcol=NULL, 
    border="black", all.labels=FALSE, ask=if(mixed & write.plot=="show") 
    TRUE else FALSE )
}
\arguments{
\item{x}{
for the \code{mvabund} method \code{x} specifies the data from
which the boxplots are to be produced. 
This can be either a numeric vector, or a single list containing such vectors. 
Additional unnamed arguments specify further data as separate vectors 
(each corresponding to a component boxplot). 
NAs are allowed in the data. \cr
For the default method, unnamed arguments are additional data vectors 
(unless \code{x} is a list when they are ignored), and named arguments are
arguments and graphical parameters to be passed to in addition to the ones given
by argument pars (and override those in pars). \cr
For the \code{mvformula} method, a formula, such as \code{y ~ grp}, where
y is a numeric mvabund object of data values to be split into groups according
to the grouping variable grp (a factor). 
}
\item{y}{for the \code{mvabund} method \code{y} can be an additional
\code{mvabund} object, if \code{x} isn't a list.}
\item{\dots}{
for the \code{mvformula} method, named arguments to be passed to
the \code{plot.mvformula} method. Some arguments that are available
for the \code{mvabund} method, are not available in \code{plot.mvformula}
and can therefore not available in the \code{mvformula} method.

For the \code{mvabund} method, unamed arguments are additional data of vectors or matrices or \code{mvabund} objects, (unless \code{x} is a list when they are ignored),and named arguments are arguments and graphical parameters to be passed in addition to the ones given by argument \code{pars} 
(and override those in \code{pars}).
}
\item{data}{for \code{boxplot.mvformula} an optional data.frame that
contains the variables of the mvformula \code{x}}
\item{subset}{for \code{boxplot.mvformula} an optional subset of the data
that should be used}
\item{range}{this determines how far the plot whiskers extend out from the box.
If range is positive, the whiskers extend to the most extreme data point which
is no more than range times the interquartile range from the box. A value of
zero causes the whiskers to extend to the data extremes.}
\item{width}{a vector giving the relative widths of the boxes making up the
plot.}
\item{notch}{if \code{notch} is \code{TRUE}, a notch is drawn in
    each side of the boxes.  If the notches of two plots do not
    overlap this is \sQuote{strong evidence} that the two medians differ
    (Chambers \emph{et al.}, 1983, p. 62).  See \code{boxplot.stats}
    for the calculations used.}
\item{outline}{if \code{outline} is not \code{TRUE}, the outliers are
    not drawn (as points whereas S+ uses lines).}
    % the argument name is most ugly but S+ compatible
\item{names}{only available for the \code{mvabund} method:
group labels which will be printed under each boxplot.}
\item{boxwex}{a scale factor to be applied to all boxes.  When there
are only a few groups, the appearance of the plot can be improved
by making the boxes narrower.}
% \item{plot}{if \code{TRUE} (the default) then a boxplot is
%     produced.  If not, the summaries which the boxplots are based on
%     are returned.}
\item{border}{an optional vector of colors for the outlines of the boxplots.
The values in \code{border} are recycled if the length of \code{border} is
less than the number of plots. A list of the possible colors can be obtained
with the function \code{colors}. }
\item{col}{if \code{col} is not \code{NULL} it is assumed to contain colors
    to be used to colour the bodies of the box plots. By default they
    are in the background colour.}
% \item{log}{character indicating if x or y or both coordinates should
%     be plotted in log scale.
% Note, that for the \code{mvabund} method \code{log} should NOT be used
% and for the \code{mvformula} method transformations should be used
% inside the formula, if wanted.
% Instead \code{transformation} can be specified. }
\item{pars}{a list of (potentially many) more graphical parameters,
    e.g., \code{boxwex} or \code{outpch}; these are passed to
    \code{bxp}; for details, see there.}
% \item{horizontal}{logical indicating if the boxplots should be
%     horizontal; default \code{FALSE} means vertical boxes.}
% \item{add}{logical, if \code{TRUE} \emph{add} boxplot to current plot.}
\item{at}{only available for the \code{mvabund} method:
numeric vector giving the locations where the boxplots should
    be drawn;
    defaults to \code{1:n} where \code{n} is the number of boxes.}
\item{xlab}{a title for the x axis, see \code{plot}.}
\item{ylab}{a title for the y axis, see \code{plot}.}
\item{main}{an optional main, a character value or a character vector (only if
two mvabund.objects are passed).}
\item{pch}{this argument is not used.}
\item{fg}{the color of the axis, default=grey.}
\item{las}{the style of axis labels.}
\item{write.plot}{
possible values: "show" (default) = plot is shown, or "eps", "postscript", "pdf", "jpeg",
"bmp", "png" to save the plot in this format, plot is then not shown.}
\item{filename}{character, pathname where the file should be stored. 
Default location is R's working directory.}
\item{n.vars}{the number of variables to include in the plot.}
\item{overall.main}{a character to display as title for every window.}
\item{var.subset}{a numeric vector of indices indicating which variables of
the mvabund.object should be included on the plot.}
\item{transformation}{
an optional transformation, (ONLY) for the \code{mvabund} method.
Note, that for the \code{mvabund} method \code{transformation} must be used 
instead of \code{log}.\cr
Available values are:\cr
"no" = untransformed, "sqrt"=square root transformed,
"log" (default)=log(Y/min+1) transformed, "sqrt4" =4th root transformed.
}
\item{scale.lab}{argument for limits and labelling of tickmarks: \cr
"r" = R's default limits are used \cr
"s" (default) = standard limits beginning at 0 are used.
}
\item{t.lab}{argument for limits and labelling of tickmarks: \cr
"o" (default) = \dQuote{original} values of y are visible in the labels, e.g.
for transformation="sqrt", then the tickmark labels are based on x (and y) not
sqrt(x) (and sqrt(y)) \cr
"t" = transformed  values of y are visible in the labels.}
\item{mfrow}{layout of the plotting window, either a pair of c(rows, columns)
per window or an integer giving the prefered number of plots per window.}
\item{mfcol}{like \code{mfrow}, but the order of drawing is by row and not by
column. If \code{mfcol} is passed, \code{mfrow} is ignored.}
% \item{checks}{logical, if applicable: whether it should be checked if
% the variables of two \code{mvabund} objects match.}
\item{ ask }{logical, logical. If \code{TRUE} the user is asked for input,
before a new figure is drawn.}
% \item{line.col}{optional, the color of the added line if applicable.}
\item{all.labels}{logical, for \code{boxplot.mvformula}, whether all labels
should be drawn or just the labels on top and on the right side of each window \cr
(this only makes a difference if there are several plots drawn per window) }
}
\details{ 
The function \code{boxplot.mvabund} allows simultaneous construction of many
variables on a single figure. Thus a good comparative overview about the
distribution of abundances for several species can be obtained. \cr
There are several ways in which this function can be used. 
If one \code{mvabund} object, either named \code{x} or \code{y} 
or not names, is passed, it will be drawn on one plot and abundances can be
compared over several variables.\cr
If two \code{mvabund} objects, named \code{x} and \code{y} are 
passed for plotting, they will be shown on
one plot, showing for each species the abundances of both objects directly
one below the other. \cr
If more than two \code{mvabund} objects are passed, each of them will be
plotted separately. \cr
Additionally, it is possible to specify \code{x} as a list of \code{mvabund} objects.
Each of them will be plotted separately and any further \code{mvabund} data will 
be ignored, regardless if it is specified as \code{y} or unnamed.

The function \code{boxplot.mvformula} can be used to draw boxplots of a \code{mvabund}
object in dependence of explanatory variables. The explanatory variables can be both 
numerical values as well as factor variables. If the formula contains both of them,
there will be separate plots for the terms with numerical values and the terms
with factor variables, displayed on separate windows.

The arguments \code{plot}, \code{varwidth} and \code{add}, which are availabe in the default method of \code{boxplot}, are not available for the \code{mvabund} and \code{mvformula} methods. The argument \code{horizontal} is not available for the \code{mvabund} method. \cr
A number of other arguments like \code{at} and \code{names} are only available for the 
\code{mvabund} method.
}
\section{Warning}{
The argument \code{log}, that is available in most plotting functions can not be used 
for plotting \code{mvabund} or \code{mvformula} objects. Instead use \code{transformation} for the \code{mvabund} method and for the \code{mvformula} method
include transformations in the formula.
}

\value{In contrast to the default method (boxplot.default) nothing will be returned.
These functions are only used for drawing the plots.
}
\references{
Warton, D. I. ( )
\emph{Raw data graphing: an informative but under-utilised tool
for the analysis of multivariate abundances}, , .
}
% \note{  }
\author{
Ulrike Naumann, Stephen Wright and David Warton <David.Warton@unsw.edu.au>.
}
\seealso{
\code{\link{plot.mvabund}}.
}
\examples{
require(graphics)

#### Basic Use ####
data(spider)
spiddat <- spider$abund
X <- spider$x

## Create the mvabund object:
spiddat <- mvabund(spiddat)

## Draw a boxplot for a mvabund object:
boxplot(spiddat)

## the same plot could be done by
plot(spiddat,type="bx")


#### Advanced Use ####
data(solberg)
solbdat <- solberg$abund
treatment<- solberg$x

# create pch type and colour vectors
treat.pch <- treat.col <- unclass(treatment)

# Boxplot for data
plot.mvabund(x=solbdat,y=treatment,type="bx",
             main="BoxPlot of The 12 Highest Abundant Species", 
             xlab="Abundance [sqrt scale]",ylab="",
             transformation="sqrt",t.lab="o",shift=TRUE)

}
\keyword{hplot}

