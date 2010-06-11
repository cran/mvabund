\name{default.plot.mvabund}
\alias{default.plot.mvabund}
\alias{default.plot.mvformula}
\title{Plot Multivariate Abundance Data and Formulae (Default Function Call)}
\description{
Produces a range of plots for visualising multivariate abundance data and 
its relationship to environmental variables, including: dot-plots and boxplots 
for different levels of a factor
stacked by response variable; comparative dot-plots and boxplots for different 
levels of a factor, stacked by response variable; scatterplots of abundances 
against a set of explanatory variables; scatterplots of pair-wise abundance 
data, e.g. from repeated measures. See details below.
}
\usage{
\method{default.plot}{mvabund}(x, y, type="p", main,xlab, ylab,col= if(type=="bx") "white" else "black", 
  fg= "grey", pch=1, las=1, write.plot="show", filename="plot.mvabund", 
  n.vars= min(12,NCOL(x)), overall.main, var.subset=NA,	subset=NA, 
  transformation="log",scale.lab="ss",t.lab="o", mfrow=if(!two.objects | (type=="bx")) { 1 } 
  else if(write.plot=="show") min(25, n.vars) else min(9,n.vars), 
  mfcol=NULL,	shift=TRUE,	border= if(two.objects & length(col)==1) { c("red","blue")} 
  else "black", checks = TRUE, add.line=FALSE, line.col="black", keep.window=FALSE, 
  ask=TRUE,\dots)

\method{default.plot}{mvformula}(x,y=NA, type="p", main="", xlab=NULL, ylab="Abundances",
  col=if(type=="bx") "white" else "black", fg= "grey", pch=1, las=1, 
  write.plot="show", filename="plot.mvabund", 
  n.vars= if(any(is.na(list(var.subset)))) 12 else length(var.subset),
  overall.main="", data = parent.frame(), var.subset=NA, 
  xvar.select=TRUE, n.xvars=if(any(is.na(xvar.subset))) 
  min(3, sum(!is.interaction)) else length(xvar.subset), 
  xvar.subset = NA, transformation.formula="log", scale.lab="ss", t.lab="t",
  mfrow=c( min(5,n.vars), min(3, n.xvars[xvar.select])), 
  mfcol=NULL, shift=TRUE, border="black", all.labels=FALSE, 
  keep.window=FALSE, ask, \dots)
}

\arguments{
\item{x}{
for the \code{mvabund} method, x is a \code{mvabund} object.\cr
For the \code{mvformula} method, x is a \code{mvformula} object,
a Model Formula to be used.
}
\item{y}{in \code{plot.mvabund} an optional second matrix with multivariate
abundance data in \code{plot.mvformula} an optional matrix of the
independent variables to explain x.}
\item{type}{what type of plot should be drawn. Useful types are "p" for
scatterplot, "bx" for boxplot and "n" for no plot. Other types, see
\code{plot} are allowed, but usually don't give a meaningful output.}
\item{main}{the title of the plot, see \code{plot}.}
\item{xlab}{a title for the x axis, see \code{plot}.}
\item{ylab}{a title for the y axis, see \code{plot}.}
\item{col}{the color of the plotting symbols. See \code{par}, as
well as the Details.}
\item{fg}{the color of the axis, default=grey.}
\item{pch}{the plotting symbol.
See \code{points} for possible values and their interpretation.
and see Details.
}
\item{las}{the style of axis labels, see \code{par}.}
\item{write.plot}{
can be "show" (default): plot is shown, or "eps", "postscript", "pdf", "jpeg",
"bmp", "png" to save the plot in this format, plot is then not shown.
}
\item{filename}{character, pathname where the file should be stored. 
Default location is R's working directory.}
\item{n.vars}{the number of variables to include in the plot.}
\item{overall.main}{a character to display as title for every window.}
\item{data}{only used for plot.mvformula: a data.frame (or list) from
which the variables in formula should be taken.}
\item{var.subset}{a numeric vector of indices indicating which variables of the
mvabund.object should be included on the plot.}
\item{subset}{a vector specifying a subset of cases, if any values in subset
are NA, the argument is discarded.}
\item{transformation}{an optional transformation, if no formula is given,
"no" = untransformed, "sqrt"=square root transformed,"log" (default)=log(Y/min+1) transformed, "sqrt4" =4th root transformed.\cr
Note that if \code{plot.mvabund} is called explicitly, and two data objects 
supplied, none of which is a \code{mvabund} object, then \code{plot.mvformula}
will be called (See Details). The argument \code{transformation} is then NOT 
available. 
}
\item{transformation.formula}{an optional transformation passed to \code{plot.mvformula}: "no" = untransformed, "sqrt"=square root transformed, "log" (default)=log(Y/min+1) transformed, "sqrt4" =4th root transformed.\cr
}
\item{xvar.select}{whether only a subset of x variables should be plotted or all.}
\item{n.xvars}{the number of the most relevant x variables that should be
plotted, is not used if \code{xvar.select = FALSE}.}
\item{xvar.subset}{
a subset of x variables that should be plotted, is not used
if \code{xvar.select = FALSE}.
}
\item{scale.lab}{argument for limits and labelling of tickmarks: \cr
"r" = R's default limits are used (only in \code{plot.mvabund}) \cr
"s" (default) = standard limits beginning at 0 are used. If two \code{mvabund}
objects are drawn, the x and y-axis have the same scale.
"ss" = (if applicable) Like "s", however the axes have the same scale
for all abundances.
}
\item{t.lab}{argument for limits and labelling of tickmarks: \cr
"o" (default) = "original" values of y are visible in the labels,  \cr
   e.g. for transformation="sqrt", then the tickmark labels are based on x
   (and y) not sqrt(x) (and sqrt(y)) \cr
"t" = transformed  values of y are visible in the labels.}
\item{mfrow}{layout of the plotting window, either a pair of c(rows, columns)
per window or an integer giving the prefered number of plots per window.}
\item{mfcol}{like \code{mfrow}, but the order of drawing is by row and not by
column, if mfcol is passed, mfrow is ignored.}
\item{shift}{logical, whether overlapping points in dotplots should be shifted,
so that every point is visible. Is not used in scatterplots.}
\item{border}{an optional value of colors for the outlines of the boxplots.}
\item{checks}{logical, if applicable: whether there should be checks if 
the variables of two \code{mvabund} objects match.}
\item{add.line}{logical, (default is \code{FALSE}), whether a 'diagonal'
line ((0,0) to (max,max)) is added to the two-mvabund-objects-scatterplot. }
\item{line.col}{optional, the color of the added line if applicable, can be
specified in the same way as color.}
\item{all.labels}{only for in \code{plot.mvformula}: logical, whether all
labels for the x- and y-axis are drawn or whether only the
x-axis labels at the bottom graphic and only the y-axis labels for the
graphic on the right side of the window are produced.}
\item{keep.window}{whether the specifications of mfrow, mfcol and windowsize
are kept, if \code{write.plot} is \code{"show"}. if \code{write.plot} is not
\code{"show"}, \code{keep.window} will not be used.
If \code{write.plot} is \code{"show"} and there is still 'space' for more plots
on the current graphic device, then any following plot(s) are drawn
on the interactive graphic device that is currently open. However, this option
should ONLY be used for adding simple plots or a legend.
It should NOT be expected to work properly for higher level plotting functions
or functions that produce multiple plots, for example another call to
\code{plot.mvabund}.}
\item{ ask }{logical, whether user input should be awaited before drawing
new windows.}
\item{\dots}{arguments to be passed to or from other methods.}
}
\details{
The function \code{plot.mvabund} produces plots for the visualisation of
multivariate abundance data and their relationships to environmental variables.
The approach taken is to separately plot the relationship between each response
variable and environmental
variables, that is, to visualise the marginal distribution, as in Warton (2008). 
Three main types of plot that can be produced:

(1) Dot-plots or boxplots stacked along the y-axis by response variable. If a
factor is given, comparative dot-plots/boxplots are produced, comparing
abundances across each factor level. This type of plot is produced when one
multivariate abundance dataset is given as an input argument, either on its own,
or together with a factor, as in the examples using the solberg dataset below.

(2) Scatterplots of multivariate abundances against environmental variables,
with separate plots for separate response variables. This type of plot is
produced when one multivariate abundance dataset is given as an input argument
together with an environmental variable or a set of environmental variables.

(3) Scatterplots of a paired sample of multivariate abundances. This type of
plot is produced when two multivariate abundance datasets are given as input
arguments, and their size and variable names match each other. It is up to the
user to ensure that the rows match for these two datasets.

There are several methods for calling \code{plot.mvabund}:  \cr
(a) plot.mvabund("matrix", ...) \cr
The multivariate abundances are stored
in the data matrix. Including an optional second argument determines
whether a plot of type (1) is produced (if no second argument or if it is a factor), 
or a plot of type (2) (if one or a set of environmental variables is given), or a 
plot of type (3) (if a second matching multivariate abundance dataset is given). \cr
Instead of a matrix, \code{mvabund} objects can be used. \cr
(b) plot("mvabund object", ...) \cr
You can define mvabund objects using the function 
\code{mvabund}. Then the behaviour of the plot function is the same as 
\code{plot.mvabund} above. \cr
(c) plot.mvformula("response"~"terms") \cr
The first of these two objects must be the multivariate abundances, which can be
either a matrix or a \code{mvabund} object. The terms determine the type of
plot produced. The terms can be either a single vector or matrix or 
a number of vectors or matrices, separated by \code{+}. 
Compare \code{formula} for further details on the specification
of the terms. \cr
(d) plot("mvformula object") \cr
You can define mvformula objects using the function \code{mvformula}.
Note that the response cannot be a data frame object.

For plots of type (3) above, you must use method (a) or (b).
Plot methods (c) and (d) require that both the response and explanatory
variables are specified, i.e. formulas like '~x' or 'y~1' cannot be plotted.

See below examples to see how each of these methods is applied.

Multivariate abundance datasets typically have many variables,
more than can be visualised in a single window, so by default plot.mvabund
subsets abundance variables (and where appropriate, environmental variables).
By default the 12 most abundant variables are plotted (determined on transformed
 variables if the response is transformed in the mvformula method),
although this setting can be controlled via the argument \code{n.vars}, and the
variables included in the subset to be plotted can be controlled via
\code{var.subset}. It is possible for example to plot the abundance variables
most significantly associated with an environmental variable,
as in the Solberg example below.

To produce boxplots rather than dot-plots in type (1) plots, use the argument
\code{type="bx"}.

For type (2) plots, if only one environmental variable is specified, plots for 
different abundance variables are arranged in a rectangular array (different 
abundance variables in different rows and columns). If however more than one
environmental variable is specified, different columns correspond to different
environmental variables (and different abundance variables in different rows).
If more than 3 environmental variables are specified, the 3 will be selected
that maximise average R^2 when \code{manylm} is applied (using the subset
selection function \code{best.r.sq}). To switch off this subset selection, set
\code{xvar.select=FALSE}, or choose your own subset of environmental variables
using \code{xvar.subset}.

To control the appearance of points on dot-plots and scatterplots, usual
arguments apply (see \code{par} for details). The plotting symbols \code{pch}
and their \code{color} can be a vector, and if the plot function is called via
a mvformula object, it can also be a list, where the list elements corresponds
to the symbols / colors used in the plots for different
independent variables. 

If some of the formula terms are factor variables, these will be drawn in
boxplots.
Note, that the plots produced by \code{plot.mvformula}
depend on whether the first independent variable is a factor or not. 
See the examples for the different possibilities of boxplots that can be
produced.

If two objects are passed and only one of them is an \code{mvabund} object,
the resulting plots will be the same as if a formula was supplied, having the
\code{mvabund} object as response variable.\cr
If both objects are not \code{mvabund} objects, it will be tried to guess which
one of them is the response. The following logic applies:
If \code{y} is not a \code{data.frame}, it will be assumed that \code{y} is the
response. Note that \code{y} is the second object, if argument names are not
supplied.
If \code{y} is a \code{data.frame} and \code{x} is not a \code{data.frame},
it will be assumed that \code{x} is the response. If both objects are
data frames an error will result, as the function is designed for \code{mvabund}
objects!

The argument \code{shift} controls whether or not points are shifted on dotplots 
so that they do not overlap. This argument is ignored for boxplots and
scatterplots (type (2) or type (3) graphs).
}
\section{Warning}{
The argument \code{log}, that is available in lots of plotting functions can not
be used for plotting \code{mvabund} or \code{mvformula} objects. Instead use
\code{transformation} for the \code{mvabund} method and for the
\code{mvformula} method include any transformations in the formula.
}
\references{
Warton D.I. (2008). Raw data graphing: an informative but under-utilized tool
for the analysis of multivariate abundances. \emph{Austral Ecology} 33(3), 290-300.
}
\author{
Ulrike Naumann, Stephen Wright and David Warton <David.Warton@unsw.edu.au>.
}
\seealso{
\code{\link{boxplot.mvabund}}, \code{\link{meanvar.plot}}, \code{\link{plot.manylm}}, 
\code{\link{plot.manyglm}}.
}
\examples{

require(graphics)

############################
## Some ``type (1)" plots ##
############################

data(solberg)
solbdat <- solberg$abund
treatment<- solberg$x

## Plot a multivariate dataset (Species vs Abundance)
plot.mvabund(solbdat)

## Alternatively, the plot command could be used directly if spiddat is 
## defined as an mvabund object:
solbmvabund <- mvabund(solbdat)
plot(solbmvabund)

## Draw an mvabund object in a boxplot, but using the 20 most abundant 
## variables in the plot, using the square root transform, and adding 
## coloured axes and title:
plot.mvabund(solbdat, n.vars=20, type="bx", transformation="sqrt", 
fg="lightblue", main="Solberg abundances", col.main="red")

## Plot Species (split by treatment) vs Abundance
plot(solbmvabund,treatment)

## This can also be produced using
plot.mvabund(solbdat, treatment)

## or
plot.mvformula(solbmvabund~treatment)


## To use plot.mvabund to plot only the variables with P-values less than 0.1:
lm.solberg <- manylm(log(solbmvabund+1)~treatment)
anova.solb <- anova(lm.solberg, p.uni="unadjusted")
pj = anova.solb$uni.p

plot(solbmvabund~treatment, var.subset=pj<0.1)

## Or to plot only the 12 most significant variables, according to their 
## univariate ANOVA P-values:
pj.sort = sort(pj, index.return=TRUE)
plot(solbmvabund~treatment, var.subset=pj.sort$ix[1:12])



############################
## Some ``type (2)" plots ##
############################
#load and convert data
data(spider)
spiddat <- mvabund(spider$abund)
spidx <- mvabund(spider$x)

#create labels vectors
pch.vec <- as.numeric(spidx[,3]<2)
pch.vec[pch.vec==0] <- 3

#Scale the soil water variable
soilWater <- spidx[,1]

#Create the Table for the main titles of each plot
title <- c("\n\nAlopecosa accentuata", "\n\nAlopecosa cuneata",
           "\n\nAlopecosa fabrilis", "\n\nArctosa lutetiana", 
           "\n\nArctosa perita", "\n\nAulonia albimana", 
           "\n\nPardosa lugubris", "\n\nPardosa monticola", 
           "\n\nPardosa nigriceps", "\n\nPardosa pullata",
           "\n\nTrochosa terricola", "\n\nZora spinimana")

#Plot Species Abundance vs Environmental variable
plot.mvformula(log(spiddat+1)~soilWater, main=title, xlab="log(\% Soil Moist)", ylab="Abund [log(y+1) scale]", overall.main="Species Abundance vs \% Soil Moisture", fg="grey", col=pch.vec, pch=pch.vec, las=1, scale.lab="ss", mfrow=c(4, 3))


#Add a Margin
par(xpd=NA)
legend("topright",pch=c(1,3),col=c(1,3),legend = c("few twigs", "many twigs"), cex = 0.8, inset=c(-0.05,-0.24))


############################
## Some ``type (3)" plots ##
############################

##Plot 1981 Abundance vs 1983 Abundance
data(tikus)
year <- tikus$x[,1]
tikusdat <- mvabund(tikus$abund)
site <- tikus$x[,2]

plot(tikusdat[year==81,], tikusdat[year==83,], col.main="blue",
  	xlab="1981 abundance", ylab="1983 abundance")


}
\keyword{hplot}

