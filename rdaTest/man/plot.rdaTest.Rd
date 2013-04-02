\encoding{UTF-8}
\name{plot.rdaTest}
\alias{plot.rdaTest}
\alias{multi.factor}
\alias{ellipse}
\alias{larger.frame}
\title{ Canonical redundancy analysis triplots }
\description{
Function \code{plot.rdaTest} produces triplots from the output list of \code{\link{rdaTest}}.
}
\usage{
\method{plot}{rdaTest} (x, xax=1, yax=2, scaling=1, graph.type="notchosen", plot.sites=TRUE, 
plot.spe=TRUE, plot.env=TRUE, binary=NULL, height=6, width=6, xlim=NULL, 
ylim=NULL, ell=NULL, alpha=0.95, ell.axis=FALSE, mul.spe=0.8, mul.env=0.90, 
mul.text=0.10, pos.ell=NULL, pos.bin=NULL, pos.sites=1, label.sites=TRUE, 
label.spe=TRUE, label.env=TRUE, label.ell=TRUE, col.env="blue", col.spe="red", 
col.site="black", col.ell="black", lty.env=1, lty.spe=1, lty.ell=1, lty.axis=2, 
cex=1, cex.lab=1, cex.axis=1, lwd=1, len=0.1, saveplot=FALSE, path=NULL, 
mar.perc=NULL, interior.mar.perc=0.05, select.spe=NULL, ...)
}

\arguments{

  \item{********** }{ General plotting parameters }

  \item{x }{ The rdaTest output object. The most simple call to this function is to provide only the output object name of function \code{\link{rdaTest}}. Ex.: plot.rdaTest(toto) or plot(toto). }
  \item{xax, yax }{ Canonical axes that will be used for the triplot. These value can be negative; if so, the plot axes will be reversed. }
  \item{scaling }{ Allows the user to choose between the two scaling types, 1 or 2, for the triplots: scaling = 1: distance triplot (default); scaling = 2: correlation triplot. }
  \item{graph.type }{ Either "F" or "Z" -- "F" = site scores, "Z" = fitted site scores. }
  \item{plot.sites }{ If TRUE, the sites will be plotted as small circles. }
  \item{plot.spe }{ If TRUE, the species (or other response variables) will be plotted. }
  \item{plot.env }{ If TRUE, the explanatory (environmental) variables will be plotted. }
  \item{binary }{ A vector specifying which of the environmental variables are binary; the corresponding data points will be averaged and represented by a symbol instead of an arrow, as in Canoco. }
  \item{saveplot }{ If \code{saveplot = FALSE}, the graph will only appear in a window (quartz). If \code{saveplot = TRUE}, no graph will be plotted on the screen; a pdf file will created instead. }
  \item{path }{ Complete path to the file in which the plot will be saved; example: \code{~/Desktop/toto.pdf} }
  \item{width, height }{ Size of the triplot in inches. }
  \item{xlim, ylim }{ Vectors describing the minimum and maximum values of the plotted region. }
  \item{lwd }{ Line width for the axes, arrows, and ellipses. }
  \item{len }{ Length of the arrow heads. }
  \item{mar.perc }{ Allows one to increase or decrease the size of the outer margin of the graph, using a fraction of the total width and height of the graph; the value must be between 0 and 1. }
  \item{interior.mar.perc }{ Allows one to increase the interior margin of the plot, i.e. the space between the outermost points and the frame, using a fraction of the total size of the frame; the value must be between 0 and 1. Default value: \code{interior.mar.perc=0.05}. }
  \item{select.spe }{ Vector containing a selection of the species, to  be drawn in the biplot. Draw all species if select.spe=NULL (default value). }
  \item{...}{ Other graphical arguments passed to the function. }

  \item{********** }{ }
  \item{********** }{ Parameters controlling the confidence ellipses around groups of points }

  \item{ell }{ Draw confidence ellipses around groups of points according to a vector defining the group assignments of the objects; see section Examples below. }
  \item{alpha }{ Confidence region of the ellipses (e.g. 0.90, 0.80, etc.) Default value: 0.95. }
  \item{lty.ell }{ Drawing type for ellipses; see \code{lty} in the help documentation file of \code{\link{par}}. }
  \item{col.ell }{ The ellipses and group identification text color; see \code{lty} in the help file of \code{\link{par}}. }
  \item{ell.axis }{ If TRUE, draw the major and minor axes of the ellipses. }
  \item{pos.ell }{ Offset the group names relative to the ellipse centers: 1=bottom, 2=left, 3=top, 4=right. }

  \item{********** }{ }
  \item{********** }{ Parameters for arrows and text }

  \item{label.sites }{ If TRUE, labels are added to the site symbols. }
  \item{label.spe }{ If TRUE, labels are added to the species arrows. }
  \item{label.env }{ If TRUE, labels are added to the environmental variable arrows. }
  \item{label.ell }{ If TRUE, labels (numbers) are added to the group ellipses. }
  \item{col.env }{ Color for the environmental variable arrows and text. }
  \item{col.spe }{ Color for the species arrows and text. }
  \item{col.site }{ Color for the site points and text. }
  \item{lty.env }{ Environmental variable arrows drawing type. See \code{lty} in the help documentation file of \code{\link{par}}. }
  \item{lty.spe }{ Species arrows drawing type.  See \code{lty} in the help documentation file of \code{\link{par}}. }
  \item{lty.axis }{ Drawing type for axes, full lines, dotted lines, etc. }

  \item{cex }{ A numerical value giving the amount by which plotting text and symbols should be magnified relative to the default. See \code{lty} in the help documentation file of \code{\link{par}}. }
  \item{cex.lab }{ Change to the font size of the site, species, and environmental variable labels, with respect to a standard of 1. }
  \item{cex.axis }{ Change to the font size of the axis labels, with respect to a standard of 1. }
  \item{mul.spe }{ Multiplier based on the range of the plot axes, which changes to the length of the species arrows proportionally; 1 = maximum relative length in the set of arrows. }
  \item{mul.env }{ Same, for the explanatory variables. }
  \item{mul.text }{ Text position at the end of each arrow plus a percentage; default is 10 percent. }
  \item{pos.bin }{ Offset the names of the binary variables: 1=bottom, 2=left, 3=top, 4=right. }
  \item{pos.sites }{ Offset the site labels: 1=bottom, 2=left, 3=top, 4=right. }
}

\details{
A triplot of the results of the canonical redundancy analysis (RDA) is drawn using either scaling 1 or scaling 2. 

Scaling type 1 triplot uses matrices Z or F for objects and U for variables; notation as in Legendre and Legendre (1998, Section 11.1.1). This projection preserves the Euclidean distances among the objects.

Scaling type 2 triplot uses matrices ZSc2 or FSc2 for objects and USc2 for variables. This projection preserves the correlations among the response variables' fitted values.
}

\author{ Sebastien Durand and Pierre Legendre, Universite de Montreal }

\references{ 
Aart, P. J. M. (van der) and N. Smeenk-Enserink. 1975. Correlations between distributions of hunting spiders (Lycosidae, Ctenidae) and environmental characteristics in a dune area. Neth. J. Zool. 25: 1-45.

Legendre, P. and Legendre, L. 2012. Numerical Ecology. 3rd English ed. Elsevier, Amsterdam. 
}

\seealso{ \code{\link{rdaTest}} }

\examples{
## Analysis of Table 11.3 of Legendre and Legendre (1998, p. 590)  
data(Table.11.3)
Y = Table.11.3[,1:6]
X = Table.11.3[,7:10]
result = rdaTest(Y, X, test.F=TRUE, nperm=999)
plot(result, graph.type="Z")

## Draw 3 species only in the biplot
vec = c(1,3,6)
plot(result, graph.type="Z", select.spe=vec)

## Reverse orientation of abscissa (parameter xax) in the plot
plot(result, xax=-1, yax=2, mul.spe=0.90, mul.env=0.70, mul.text=0.10, 
scaling=1, graph.type="F", mai.perc=0.15, pos.sites=3) 

## Draw confidence ellipses around groups of points, spider data 
## (Aart and Smeenk-Enserink, 1975)
data(spiders.spe)
data(spiders.env) 
# Hellinger transformation of the species abundance data prior to RDA
library(vegan)
spiders.hel = decostand(spiders.spe, "hellinger")
result.spiders = rdaTest(spiders.hel, spiders.env, test.F=TRUE, nperm=999)
# Vector vec divides the sites into three a priori groups
vec = c(2,1,2,1,1,3,1,3,2,2,2,2,1,3,3,3,3,3,3,3,3,2,2,2,3,2,2,2)
plot(result.spiders, graph.type="Z", ell=vec, lty.ell=3, pos.sites=4, 
mar.perc=0.15)

## Plot the sites and environmental variables, but not the species
plot(result.spiders, graph.type="Z", plot.spe=FALSE, ell=vec, lty.ell=3, 
pos.sites=4, mar.perc=0.15)
}

\keyword{ multivariate }