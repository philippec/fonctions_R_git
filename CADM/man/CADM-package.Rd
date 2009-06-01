\name{CADM-package}
\alias{CADM-package}
\docType{package}
\title{ Congruence among distance matrices }
\description{
Compute and test the coefficient of concordance among several distance matrices. Carry out a posteriori tests on the contribution of individual matrices to the overall concordance of the group.
}
\details{
\tabular{ll}{
Package: \tab CADM\cr
Type: \tab Package\cr
Version: \tab 1.1\cr
Date: \tab 2009-05-19\cr
License: \tab Copyleft\cr
LazyLoad: \tab yes\cr
}
Function \code{\link{CADM.global}} compute and test the coefficient of concordance among several distance matrices through a permutation test.

Function \code{\link{CADM.post}} carries out a posteriori tests of the contributions of individual distance matrices to the overall concordance of the group.

Use in phylogenetic analysis: to identify congruence among distance matrices (D) representing different genes or different types of data. Congruent D matrices correspond to data tables that can be used together in a combined phylogenetic or other type of multivariate analysis.
}
\author{
Pierre Legendre, Universite de Montreal

Maintainer: Pierre Legendre <pierre.legendre@umontreal.ca>
}
\references{
Campbell, V., P. Legendre and F.-J. Lapointe. 2009. Assessing congruence among ultrametric distance matrices. Journal of Classification (In press).

Campbell, V., P. Legendre and F.-J. Lapointe. Performance of the congruence test among distance matrices in phylogenetic analysis. (Submitted MS).

Legendre, P. and F.-J. Lapointe. 2004. Assessing congruence among distance matrices: single malt Scotch whiskies revisited. Australian and New Zealand Journal of Statistics 46: 615-629. 

Legendre, P. et F.-J. Lapointe. 2005. Congruence entre matrices de distance. P. 178-181 in: Makarenkov, V., G. Cucumel et F.-J. Lapointe [eds] Comptes rendus des 12emes Rencontres de la Societe Francophone de Classification, Montreal, 30 mai - 1er juin 2005.
}
\keyword{ package }
\seealso{ \code{\link{kendall.W}}, \code{\link[ape:ape-package]{ape}} }
\examples{
data(mat3)
res.global <- CADM.global(mat3, 3, 9, nperm=999)
res.post   <- CADM.post(mat3, 3, 9, nperm=999, mantel=TRUE)
}
