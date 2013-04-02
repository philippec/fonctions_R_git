\name{test.axes}
\alias{test.axes}
\title{Test individual axes in RDA}

\description{
Test of significance of the individual canonical axes in simple or partial RDA using the marginal and forward methods, in the presence of a matrix of covariables W. 
}

\usage{
test.axes(Y, X, W=NULL, scale.Y=FALSE, nperm=499, seed=NULL, verbose=FALSE)
}

\arguments{
  \item{Y }{ An (n x p) response (e.g. site-by-species) data table. }

  \item{X }{ An (n x m) table of explanatory variables. The number of variables (m) is recomputed after eliminating collinear explanatory variables, if any. X will be standardized by the function. }

  \item{W }{ An optional (n x q) table of covariables. If present, W will be standardized by the function. }

  \item{scale.Y }{ A logical value (TRUE or FALSE) defining if Y should be standardized (TRUE), or only centred on its column means (FALSE). }

  \item{nperm }{ Number of permutations for the tests of the canonical axes. Values of 499 or 999 are recommended for routine work. } 

  \item{seed }{ This option allows users to set the seed for the random number generator for the test of significance of the F-statistic. The seed is used by the permutation function \code{\link{sample}}. If \code{seed=NULL} (default), a random integer is drawn as the seed for the run. It will be reset to that value before the test of each canonical axis. All axes are thus tested using the same set of permutations. A fixed value, e.g. seed=12345, can be givenn by the user to compare the results of this function with that of other functions where the seed can also be set at run time. }
  
  \item{verbose }{ If \code{verbose = TRUE}, intermediate F and F.perm results are printed during testing. By default, \code{verbose = FALSE}. }
}

\details{

The 'marginal' method is the one implemented in function permutest.cca() of the \code{vegan} package. This is the function performing the tests of significance when users call the anova.cca() function after canonical analysis by the functions rda() or cca(). The 'forward' method is implemented in the program Canoco (ter Braak and Šmilauer 2002).

The function uses permutation of residuals of the reduced model (Freedman & Lane, 1983).

This function was originally written as an explicit description of the 'marginal' and 'foward' testing methods in RDA, in the presence of a matrix of covariables W. It was published in Supplement 2 of the Legendre et al. (2011) paper, where its name is test.axes.cov.R.
}

\value{ 

The function returns a list containing the following results and matrices:

\item{test.axes }{ A table containing the results of the tests of significance of the individual canonical axes, using both the marginal and forward methods. The columns are: Eigenvalue, F.marginal = F-statistic for the marginal test, F.forward = F-statistic for the forward test, P.marginal = P-value for the marginal test, P.forward = P-value for the forward test. }

\item{axes }{ A table containing the canonical axes. }
}

\author{ Pierre Legendre, Cajo J. F. ter Braak and Jari Oksanen }

\references{ 
Freedman, D., and D. Lane. 1983. A nonstochastic interpretation of reported significance levels. Journal of Business & Economic Statistics 1: 292-298.

Legendre, P. and Legendre, L. 2012. Numerical Ecology. 3rd English ed. Elsevier, Amsterdam. 

Legendre, P., J. Oksanen, and C. J. F. ter Braak. 2011. Testing the significance of canonical axes in redundancy analysis. Methods in Ecology and Evolution 2: 269-277.

ter Braak, C. J. F., and P. Smilauer. 2002. CANOCO reference manual and CanoDraw for Windows user’s Guide – Software for canonical community ordination (version 4.5). Microcomputer Power, Ithaca, New York. 
}

\seealso{ \code{\link{rdaTest}} }

\examples{
## Analysis of Table 11.3 of Legendre and Legendre (2012, p. 644)  
data(Table.11.3)
Y = Table.11.3[,1:6]
X = Table.11.3[,7:10]
out = test.axes(Y, X, nperm=999)
out$test.axes

## Analysis of the spider data  
data(spiders.spe)
data(spiders.env) 
# Hellinger transformation of the species abundance data prior to RDA
library(vegan)
spiders.hel = decostand(spiders.spe, "hellinger")
out = test.axes(spiders.hel, spiders.env, nperm=499)
out$test.axes

## Partial RDA of the spider data  
out = test.axes(spiders.hel, spiders.env[,-1], spiders.env[,1], nperm=499)
out$test.axes
}

\keyword{ multivariate }
