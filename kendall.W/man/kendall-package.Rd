\name{kendall-package}
\alias{kendall-package}
\alias{kendall.W}
\docType{package}
\title{
Kendall concordance and a posteriori tests
}
\description{
Compute and test the coefficient of concordance among several judges (variables, species). Carry out a posteriori tests on the contribution of individual variables (species) to the overall concordance of their group.
}
\details{
\tabular{ll}{
Package: \tab kendall\cr
Type: \tab Package\cr
Version: \tab 2.8\cr
Date: \tab 2009-05-19\cr
License: \tab Copyleft\cr
LazyLoad: \tab yes\cr
}
Function \code{\link{kendall.global}} compute and test the coefficient of concordance among several variables (e.g. species) through a permutation test. If several groups of variables (species) are identified a priori in the data table, coefficients of concordance will be computed for each group separately.

Function \code{\link{kendall.post}} carries out a posteriori tests of the contributions of individual variables (species) to the overall concordance of their group. If several groups of variables (species) are identified a priori in the data table, the a posteriori tests will be computed for each group separately.

Used in ecology: to discover species associations.
}
\author{
F. Guillaume Blanchet and Pierre Legendre

Maintainer: F. Guillaume Blanchet <gblanche@ualberta.ca>
}
\references{
Legendre, P. 2005. Species associations: the Kendall coefficient of concordance revisited. Journal of Agricultural, Biological, and Environmental Statistics 10: 226-245.
}
\keyword{ package }
\seealso{
\code{\link[labdsv:labdsv-package]{labdsv}}, \code{\link[vegan:vegan-package]{vegan}}
}
\examples{
library(vegan)
data(mite)
mite.hel <- decostand(mite, "hel")

# Reproduce the results shown in Table 2 of Legendre (2005), a single group
mite.small <- mite.hel[c(4,9,14,22,31,34,45,53,61,69),c(13:15,23)]
out1 <- kendall.global(mite.small)
out2 <- kendall.post(mite.small, mult="holm")

# Reproduce the results shown in Tables 3 and 4 of Legendre (2005), 2 groups
group <-c(1,1,2,1,1,1,1,1,2,1,1,1,1,1,1,2,1,2,1,1,1,1,2,1,2,1,1,1,1,1,2,2,2,2,2)
out3 <- kendall.global(mite.hel, group=group)
out4 <- kendall.post(mite.hel, group=group, mult="holm", nperm=99)
}
