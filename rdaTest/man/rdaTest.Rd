\name{rdaTest}
\alias{rdaTest}
\alias{probFrda}
\alias{FractionBySpecies}
\alias{vif}
\alias{rdaTest}
\alias{spiders.spe}
\alias{spiders.env}
\alias{Table.11.3}
\title{ Canonical redundancy analysis (RDA) }
\description{
Function \code{rdaTest} computes a simple or partial canonical redundancy analysis (RDA), with permutation tests, following the algorithm described in Legendre and Legendre (1998, Chapter 11).
}
\usage{
rdaTest(Y, X, W=NULL, scale.Y=FALSE, test.F=NULL, nperm=NULL, silent=FALSE)
}

\arguments{
  \item{Y }{ An (nxp) response (e.g. site-by-species) data table. }

  \item{X }{ An (nxm) table of explanatory variables. The number of variables (m) is recomputed after eliminating collinear explanatory variables, if any. }

  \item{W }{ An optional (nxq) table of covariables. }

  \item{scale.Y }{ A logical value (TRUE or FALSE) defining if Y should be standardized (TRUE), or only centred on its column means (FALSE).  }

  \item{test.F }{ If NULL (default), the program will ask if the user wishes to test the F-statistic. If test.F is TRUE or FALSE, no question will be asked; the program will perform the test, or not, in accordance with that indication. }

  \item{nperm}{ Number of permutations for the F-test. If \code{nperm = NULL} (default), a question will be asked by the program. If \code{nperm = 0}, permutation tests will not be computed. For large data files, the permutation test is slow since the permutation procedure is not compiled. } 

  \item{silent}{ If \code{silent = TRUE}, the basic rdaTest output is printed to the R window. See the paragraph Value. Choose \code{silent = FALSE} when the function is called by a numerical simulation function. }
}

\details{

Function \code{rdaTest} computes [partial] canonical redundancy analysis (RDA) of a data table and produces scaling 1 and scaling 2 biplots. The response variables must be standardized (scale.Y=TRUE) if they are not all expressed in the same physical dimensions. The default is \code{scale.Y = FALSE}.

Scaling type 1 biplot uses matrices Z or F for objects and U for variables; notation as in Legendre and Legendre (1998, Section 9.1). This projection preserves the Euclidean distances among the objects.

Scaling type 2 biplot uses matrices ZSc2 or FSc2 for objects and USc2 for variables. This projection preserves the correlations among the variables.

This function runs correctly with a single response variable in Y or a single explanatory variable in X. Since there is a single canonical axis in both of these cases, the first PCA axis of the residuals is added to the tables in the output list to allow the drawing of biplots.

The function prints the following information to the R window unless \code{silent = TRUE}:
	\enumerate{
    \item
	The variance inflation factors (VIF) for the explanatory variables. The value is 0 for entirely collinear variables. The covariables are not included in this calculation.
    \item
	The bimultivariate redundancy statistic (canonical R-square), as well as the adjusted R-square when there are no covariables in the analysis.
    \item
	Test of significance of the canonical relationship: the F-statistic and permutational probability, if the test of significance has been carried out. 
    \item
	The number of objects, number of response variables, and number of explanatory variables after removing collinear variables; the number of canonical eigenvalues larger than 0. 
    \item
	The total variance in matrix Y, i.e., the total sum-of-squares divided by (n-1). 
    \item
	The eigenvalues, relative eigenvalues, and the cumulative percentage of the variance of species data accounted for by the sucessive canonical eigenvalues. 
	}
}

\value{ 

The function returns a list containing the following results and matrices:

\item{VIF }{ Variance inflation factors for the explanatory variables X; the value is 0 for entirely collinear variables. The covariables are not included in this calculation. }
  \item{canEigval }{ Canonical eigenvalues. }
  \item{U }{ Canonical eigenvectors normalized to 1 (scaling 1). }
  \item{USc2 }{ Canonical eigenvectors normalized to sqrt(eigenvalue) (scaling 2). }
  \item{F.mat }{ Matrix of object scores (scaling 1). }
  \item{Z.mat }{ Matrix of fitted object scores (scaling 1). }
  \item{FSc2 }{ Matrix of object scores (scaling 2). }
  \item{ZSc2 }{ Matrix of fitted object scores (scaling 2). }
  \item{biplotScores1 }{ Biplot scores of explanatory variables (scaling 1). }
  \item{biplotScores2 }{ Biplot scores of explanatory variables (scaling 2). }
  \item{FitSpe }{ Table of cumulative fit per species (in percent) as fraction of variance of species. }
  \item{VarExpl }{ Vector of total percent fit per species after all canonical axes. }
  \item{ProbFrda }{ Probability associated with F test of the canonical relationship. }
  \item{X.mat }{ Original X matrix (required by the plotting function). }
  \item{Rsq }{ Canonical R-square. }
  }

\author{ Pierre Legendre, Universite de Montreal }

\references{ 
Aart, P. J. M. (van der) and N. Smeenk-Enserink. 1975. Correlations between distributions of hunting spiders (Lycosidae, Ctenidae) and environmental characteristics in a dune area. Neth. J. Zool. 25: 1-45.

Legendre, P. and Legendre, L. 1998. Numerical Ecology. 2nd English ed. Elsevier, Amsterdam. 
}

\seealso{ \code{\link{plot.rdaTest}} }

\examples{
## Analysis of Table 11.3 of Legendre and Legendre (1998, p. 590)  
data(Table.11.3)
Y = Table.11.3[,1:6]
X = Table.11.3[,7:10]
result = rdaTest(Y, X, test.F=TRUE, nperm=999)
plot(result, graph.type="Z")

## Run rdaTest in such a way that it prints no output to the R console
resuls <- rdaTest(Y, X, test.F=TRUE, nperm=999, silent=TRUE) 

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

## See graphical options for triplots in \code{\link{plot.rdaTest}}.
}

\keyword{ multivariate }
