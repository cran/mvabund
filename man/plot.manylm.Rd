\name{plot.manylm}
\alias{plot.manylm}
\alias{plot.manyglm}
\title{Plot Diagnostics for a manylm or a manyglm Object}
\description{
  Four plots (selectable by \code{which}) are currently available: a plot
  of residuals against fitted values, a Normal Q-Q plot, 
  a Scale-Location plot of \eqn{\sqrt{| residuals |}} against fitted values, 
  a plot of Cook's distances versus row labels.  
  By default, all of them are provided. 

  The function is not yet available for manyglm object
}
\usage{
\method{plot}{manylm}(
  x, which=1:4, caption=c("Residuals vs Fitted", "Normal Q-Q", 
  "Scale-Location", "Cook's distance"), overlay=TRUE, 
  n.vars=12, var.subset=NULL, panel=if (add.smooth) 
  panel.smooth else points, sub.caption=NULL, main="", ask, 
  ..., col=palette(rainbow(n.vars+1))[2:(n.vars+1)], 
  id.n=if (overlay) 0 else 3, 
  labels.id=rownames(as.matrix(residuals(x))),
  cex.id=0.75, qqline=TRUE, cook.levels=c(0.5, 1),
  add.smooth= if(!is.null(getOption("add.smooth")))
  {getOption("add.smooth")} else TRUE,  label.pos=c(4, 2), 
  cex.caption=1, asp=1, legend.pos= if(length(col)==1) 
  "none" else "nextplot",  mfrow= if (overlay)
  { length(which)+(legend.pos=="nextplot")} else if 
  (write.plot=="show") c(min(n.vars,3),length(which)) 
  else length(which), mfcol=NULL, write.plot="show", 
  filename="plot.mvabund", keep.window=if(is.null(c(mfrow,mfcol)))
  TRUE else FALSE, studentized= TRUE)

\method{plot}{manyglm}(
  x, which=1:4, caption=c("Residuals vs Fitted", "Normal Q-Q", 
  "Scale-Location", "Cook's distance"), overlay=TRUE, 
  n.vars=12, var.subset=NULL, panel=if (add.smooth) 
  panel.smooth else points, sub.caption=NULL, main="", ask, ...,
  col=palette(rainbow(n.vars+1))[2:(n.vars+1)],
  id.n= if (overlay) 0 else 3, 
  labels.id=rownames(as.matrix(residuals(x))),
  cex.id=0.75, qqline=TRUE, cook.levels=c(0.5, 1), 
  add.smooth= if( !is.null(getOption("add.smooth")) )
  {getOption("add.smooth")}
  else TRUE, label.pos = c(4, 2), cex.caption = 1, asp = 1,
  legend.pos= if(length(col)==1) "none" else "nextplot",
  mfrow= if(overlay) { length(which)+(legend.pos=="nextplot")}
  else if(write.plot=="show") c(min(n.vars,3) ,length(which)) 
  else length(which), mfcol=NULL, write.plot="show", 
  filename="plot.mvabund", keep.window=if(is.null(c(mfrow,mfcol))) 
  TRUE else FALSE,  studentized=TRUE)

}
\arguments{
 \item{x}{ \code{manylm} object or \code{manyglm} object, typically the 
  result of a call to \code{\link{manylm}} or \code{manyglm}.}
  \item{which}{if a subset of the plots is required, specify a subset of
    the numbers \code{1:4}.}
  \item{caption}{captions to appear above the plots}
  \item{overlay}{logical, whether or not the different variables should be
    overlaid on a single plot.}
  \item{n.vars}{the number of variables to include in the plot.}
  \item{var.subset}{the variables to include in the plot.}
  \item{panel}{panel function. The useful alternative to
    \code{points}, \code{panel.smooth} can be chosen
    by \code{add.smooth = TRUE}.}
  \item{sub.caption}{common title---above figures if there are multiple;
    used as \code{sub} (s.\code{title}) otherwise.  If \code{NULL},
    as by default, a possible shortened version of
    \code{deparse(x$call)} is used.}
  \item{main}{title to each plot---in addition to the above
    \code{caption}.}
  \item{ask}{logical; if \code{TRUE}, the user is \emph{ask}ed before
    each plot, see \code{par(ask=.)}.}
  \item{\dots}{other parameters to be passed through to plotting
    functions.}
  \item{col}{a character or integer value or vector of plotting colors, 
  different plotting colors can be used for every variable.}
  \item{id.n}{number of points to be labelled in each plot, starting
    with the most extreme.}
  \item{labels.id}{vector of labels, from which the labels for extreme
    points will be chosen.  \code{NULL} uses observation numbers.}
  \item{cex.id}{magnification of point labels.}
  \item{qqline}{logical indicating if a \code{qqline()} should be
  added to the normal Q-Q plot.}
  \item{cook.levels}{levels of Cook's distance at which to draw contours.}
  \item{add.smooth}{logical indicating if a smoother should be added to
    most plots; see also \code{panel} above.}
  \item{label.pos}{positioning of labels, for the left half and right
    half of the graph respectively, for plots 1-3.}
  \item{cex.caption}{controls the size of \code{caption}.}
  \item{asp}{the y/x aspect ratio. If the y/x aspect ratio should not be used 
   for scaling the plot, use \code{asp = NA}. 
   For the Scale-Location plot (\code{which=3}) and the plot of 
   Cooks Distance (\code{which=4}), \code{asp} is not used!}
  \item{legend.pos}{if applicable (see overlay) the position of the legend, 
   this can either be specified as a pair of coordinates or as a keyword
   specifying the position, e.g. "center", or as keywords "none" or "nextplot".}
  \item{mfrow}{a pair containing either the number of rows and columns per window 
   or the prefered number of plots per window, the order of drawing is by column;
   if \code{mfcol} is passed, \code{mfrow} is ignored. \cr
   If mfrow = NULL & mfcol = NULL, the plots will be drawn on the current window
   if one is open.}
  \item{mfcol}{like MFROW, but the order of drawing is by row.}
  \item{write.plot}{
  can be "show" (default): plot is shown, or "eps", "postscript", "pdf", "jpeg",
  "bmp", "png" to save the plot in this format, plot is then not shown.
  }
 \item{filename}{
character, pathname where the file should be stored. 
Default location is R's working directory.
}
  \item{keep.window}{logical, if TRUE and \code{write.plot=="show"}, the specifications of mfrow, 
   mfcol, oma and the size of the winodow are kept,
   if FALSE the old values will be restored, which usually means that any plot following this 
   function will be drawn on a new window.}
  \item{studentized}{logical indicating whether studentized or standardized residuals should 
   be used for plot 2 and 3.}

}
\details{
  \code{plot.manylm} is used to check the linear model assumptions that are made
  when fitting a model via \code{manylm}. Similarly, \code{plot.manyglm} checks
  the generalised linear model assumptions made when using \code{manyglm}.
  You should check the residual vs fits plot for no pattern (hence no suggestion
  of failure of key linearity and mean-variance assumptions). For manylm fits of
  small datasets, it is desirable that residuals on the normal Q-Q plot be close
  to a straight line, although in practice the most important thing is to make
  sure there are no big outliers and no suggestion of strong skew in the data.

  Some technical details on usage of this function:  \cr
   \code{sub.caption} - by default the function call - is shown as
  a subtitle (under the x-axis title) on each plot when plots are on
  separate pages, or as a subtitle in the outer margin (if any) when
  there are multiple plots per page.  % check that
  
  The \sQuote{Scale-Location} plot, also called \sQuote{Spread-Location} or
  \sQuote{S-L} plot, takes the square root of the absolute residuals in
  order to diminish skewness (\eqn{\sqrt{| E |}} is much less skewed
  than \eqn{| E |} for Gaussian zero-mean \eqn{E}).
  
  If \code{studentized=FALSE} the \sQuote{S-L}, the Q-Q, and the Residual-Leverage 
  plot, use \emph{standardized} residuals which have identical variance 
  (under the hypothesis) otherwise \emph{studentized} residuals are used.  

 If \code{mfcol} is passed, \code{mfrow} is ignored. If \code{mfrow=NULL} 
 and \code{mfcol=NULL}, par("mfrow") is used.

 Unlike other plotting functions \code{plot.manylm} and \code{plot.manyglm} 
 respectively do not have a subset argument, the subset needs to be specified 
 in the \code{manylm} or respectively \code{manyglm} function.

 For all arguments that are formally located after the position of \code{\dots},
 positional matching does not work.

 For restrictions on \code{filename} see R's help on eps/pdf/jpeg.
 Note that \code{keep.window} will be ignored if \code{write.plot} is 
 not \code{show}.
}
\author{
Ulrike Naumann and David Warton <David.Warton@unsw.edu.au>.
}
\seealso{
\code{\link{manylm}}
}
\examples{
require(graphics)

data(spider)
spiddat <- mvabund(spider$abund)
X <- spider$x

## plot the diagnostics for the linear fit of the spider data
spidlm <- manylm(spiddat~X)
plot(spidlm,which=1:2,col.main="red",cex=3,overlay=FALSE)

plot(spidlm,which=1:4,col.main="red",cex=3,overlay=TRUE)

}
\keyword{hplot}
\keyword{regression}
\keyword{multivariate}
