\name{cross}
\alias{cross}
\title{ Number of groups for hierarchical clustering by AIC and cross-validation }

\description{
  Function \code{cross} helps determine the 'best' number of clusters for a hierarchical clustering result, computed with or without a constraint of spatial or temporal contiguty. This is done by computing AICc, the Calinski-Harabasz criterion (C-H), and cross-validation.
  
}

\usage{
cross(Y, res, k1=2, k2=NULL, xv="min", xvmult=100)
}

\arguments{
  \item{Y}{ A matrix or data frame containing the input data file used in \code{\link{constrained.clust}} or \code{\link{hclust}}. }
  
  \item{res}{ The output object from \code{\link{constrained.clust}} or \code{\link{hclust}}. }

  \item{k1}{ Lowest number of clusters in the cluster membership table (default: k1=2). }

  \item{k2}{ Highest number of clusters in the cluster membership table. If k2=NULL (default), k2=(n-1)/2. This option is interesting to investigate classifications with small groups, which are found when k2 is high. Cross-validation for high values of k2 may, however, take a long time.}

  \item{xv}{ Selection of tree by cross-validation; "1se" (beware: one-s-e): gives best tree within one SE of the overall best; "min": the best tree (default); "pick": pick the tree size interactively; "none": no cross-validation. }

  \item{xvmult}{ Number of multiple cross-validation iterations. Default: xvmult=100. }

} 

\details{ The function constructs a file with the membership of all levels obtained from hierarchical clustering, and computes various statistics that may help determine the best number of groups: AICc, the Calinski-Harabasz criterion (C-H), and cross-validation residual error.

If more than 15 groups (k2 > 15) are to be tested, the user has to confirm that choice. There is also a warning stating that the cross-validation may take a long time. }

\value{ 
  \item{AIC }{ A table with 6 columns showing (1) the R-square, which is a simple information and not a criterion to identify the best partition; (2) the AICc value (the lowest value is best); (3) the Calinski-Harabasz criterion (C-H) for partitions of the objects into k1 to k2 groups (the highest value is best); (4) the P-value associated with C-H (this p-value is used here as a model criterion, not as the result of a test of hypothesis; the lowest value is best); (5) the cross-validation residual error (cvre; the lowest value is best) and (6) the  number of groups (ngr.cross) corresponding to the minimum of criterion 'cvre'. The \code{cvre} cross-validation criterion is more conservative than \code{AIC} and \code{P(C-H)}; it may thus be seen as a better criterion in most cases. }

  \item{ngr.AIC }{ Number of groups in the partition that has the lowest value of AICc. }
  
  \item{classif.AIC }{ Vector (length n) of assignment of the objects to the groups of that partition. }

  \item{ngr.P }{ Number of groups in the partition that has the lowest P-value for C-H. }
  
  \item{classif.P }{ Vector (length n) of assignment of the objects to the groups of that partition. }

  \item{cvre.and.ngr }{ Value of the cross-validation residual error (cvre) and number of groups in the partition corresponding to the pruned tree. }
  
  \item{classif.cvre }{ Vector (length n) of assignment of the objects to the groups of that partition. }

  \item{out }{ Table of assignment of the n objects to the groups in the solutions from k1 to k2 groups. }
}

\author{ Pierre Legendre, Universite de Montreal}

\references{

  Legendre, P. and L. Legendre. 1998. Numerical ecology, 2nd English
  edition. Elsevier Science BV, Amsterdam.
}

\examples{ 
###
### Artificial map data from Legendre & Legendre (1998, Fig. 13.24): n = 16

# require(vegan)
# require(spdep)
dat <- c(41,42,25,38,50,30,41,43,43,41,30,50,38,25,42,41)
coord.dat <- matrix(c(1,3,5,7,2,4,6,8,1,3,5,7,2,4,6,8, 
4.4,4.4,4.4,4.4,3.3,3.3,3.3,3.3,2.2,2.2,2.2,2.2,1.1,1.1,1.1,1.1),16,2)

# 1. Compute a Delaunay triangulation using tri2nb() of spdep 
listW <- nb2listw(tri2nb(coord.dat), style="B")
# Construct matrix links.mat 
links.mat.dat <- listw2mat(listW)

# 2. Constrained clustering using default Ward method
D.dat <- dist(dat)
res.dat <- constrained.clust(D.dat, links.mat.dat)

# 3. Cross-validation
cross.dat <- cross(dat, res.dat)
summary(cross.dat)

###
### Mite data available in "vegan": 70 sites, 35 species

# require(vegan)
data(mite)
data(mite.xy)
mite.hel <- decostand(mite, "hellinger")

# 1. Compute a Delaunay triangulation using tri2nb() of spdep 
listW <- nb2listw(tri2nb(mite.xy), style="B")
# Construct matrix links.mat 
links.mat.mite <- listw2mat(listW)

# 2. Constrained clustering using default Ward method
D.mite <- dist(mite.hel)
res.mite <- constrained.clust(D.mite, links.mat.mite)

# 3. Cross-validation (long procedure)
cross.mite <- cross(mite.hel, res.mite, k2=15)
summary(cross.mite)
cross.mite$AIC

###
### Example of spatial transect or time series

# Doubs fish data, available in ade4.
# 30 sites along the Doubs river, which runs near the France-Switzerland 
# border in the Jura mountains; 27 freshwater fish species.
# require(ade4)
data(doubs)
fish <- doubs$poi

# 1. Instead of using the coordinates of the survey sites, pretend
# that the river is a line with 30 equispaced points.
# Construct a matrix of connexions between the sampling points
# using the rank order of the points along the river.
from.to <- matrix(c(1:29, 2:30), 29, 2)
colnames(from.to) <- c("From", "To")
# ... and from it, construct the links.mat matrix
links.mat.fish <- contiguity.mat(from.to, 30)

# 2. Constrained clustering using default Ward method
fish.hel <- decostand(fish, "hellinger")
D.fish <- dist(fish.hel)
res.fish <- constrained.clust(D.fish, links.mat.fish)

# 3. Cross-validation
cross.fish <- cross(fish.hel, res.fish)
summary(cross.fish)
cross.fish$AIC
}

\keyword{ multivariate }