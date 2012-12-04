\name{constrained.clust}
\alias{constrained.clust}
\alias{const.clust}
\alias{map.hclust}
\title{ Space- and time-constrained clustering }

\description{
  Function \code{constrained.clust} carries out space-constrained
  or time-constrained agglomerative clustering from a multivariate 
  dissimilarity matrix.
}

\usage{
constrained.clust(D, links.mat=NULL, method="ward.D2", beta=-0.25, 
verbose=FALSE, target=1)
map.hclust(res, coord, k1=2, k2, title=NULL, cex=2, pos=3, col="red", ...)
}

\arguments{
  \item{D}{ A dissimilarity (distance) matrix with class \code{dist}. }
  
  \item{links.mat}{ A square matrix with 1 for connected points and 0 elsewhere. 
  The diagonal contains 0. }

  \item{method}{ The agglomeration method to be used. This should be (an 
  unambiguous abbreviation of) one of \code{"ward.D"}, \code{"ward.D2"}, 
  \code{"single"}, \code{"complete"}, \code{"average"} (= UPGMA), 
  \code{"mcquitty"} (= WPGMA), \code{"centroid"} (= UPGMC), 
  \code{"median"} (= WPGMC), or \code{"flexible"}. Method \code{"ward.D2"} 
  (default) implements the Ward (1963) clustering criterion, method \code{"ward.D"} does not (Murtagh and Legendre, submitted). }

  \item{beta}{ The beta parameter of beta-flexible clustering 
  (default: \code{beta = -0.25}). }

  \item{verbose}{ Print additional results about the relationship of a target 
  object to the other objects in the study. This option should only be used 
  when the number of objects is small. }

  \item{target}{ The number identifying an object of interest. When \code{verbose = TRUE}, additional
  information will be printed about the relationships of this object to the 
  other objects in the study during the agglomeration. }

  \item{res}{ The output object of \code{\link{constrained.clust}}, which serves as input to \code{\link{map.hclust}}. }

  \item{coord}{ File of coordinates to be used to plot the points on the maps. }

  \item{k1}{ Lowest number of clusters (default: k1=2) for maps. }

  \item{k2}{ Highest number of clusters for maps. If \code{k2=NULL} (default) or k2 > (n-1), k2 = trunc(n/2). }

  \item{title}{ Title of the pdf file, between quotes. If \code{title=NULL} (default), the pdf file is called "Constrained clustering maps.pdf" if the clustering was produced by \code{\link{constrained.clust}}, or "Unconstrained clustering maps.pdf" if it was produced by \code{\link{hclust}}.  }

  \item{cex}{ Change to the font size of the site identification numbers on the maps, with respect to a standard of 1. }

  \item{pos}{ Offset the site identification numbers on the maps: 1=underneath, 2=left, 3=above, 4=right.. }

  \item{col}{ Color of the cluster identification numbers on the maps (default: col="red"). }

  \item{...}{ Other graphical arguments passed to the \code{\link{plot}} function. }
} 

\details{ Agglomerative clustering is carried out with a constraint of spatial or temporal contiguity. This means that only the objects that are linked in matrix \code{links.mat} are considered to be candidates for clustering: the next pair of objects to cluster will be the pair that has the lowest dissimilarity value among the pairs that are linked in \code{links.mat}. 

The same rule applies during the following clustering steps, which involve groups of objects. Matrix \code{links.mat} is updated after each agglomeration step: all objects that are neighbours of each of the components that have fused are now neighbours of the newly formed cluster. 
}

\value{ 
An object of class hclust which describes the tree produced by the clustering process. The object is a list with components:

  \item{merge }{An (n-1) by 2 matrix. Row i of merge describes the merging of clusters at step i of the clustering. If an element j in the row is negative, then observation -j was merged at this stage. If j is positive then the merge was with the cluster formed at the (earlier) stage j of the algorithm. Thus negative entries in merge indicate agglomerations of singletons, and positive entries indicate agglomerations of non-singletons. }
  
  \item{height }{A set of (n-1) non-decreasing real values. The clustering height is the value of the criterion associated with the clustering method for the particular agglomeration. }
  
  \item{order }{A vector giving the permutation of the original observations suitable for plotting, in the sense that a cluster plot using this ordering and matrix merge will not have crossings of the branches. }
  
  \item{label }{Labels for each of the objects being clustered. }
  
  \item{method }{The agglomerative clustering method that has been used. }

  \item{call }{The call which produced the result. }

  \item{dist.method }{The distance that has been used to create D (only returned if the distance object has a "method" attribute). }

  \item{hcl }{An additional object, not found in the output list of \code{\link{hclust}}, which contains additional information about the constrained clustering steps: n = number of objects; len = number of distances, n*(n-1)/2; method = order number of the clustering method in the list of argument \code{method}; ia = vector of first objects in pairs; ib = vector of second objects in pairs; crit = vector of fusion level for each object at the end of the agglomeration; nn = vector (length n) with current nearest neighbor of each object, to its right in the upper-triangular matrix of distances; disnn = distance of each object to its nearest neighbor at the end of the agglomeration.}

The output of the \code{\link{map.hclust}} function is a pdf file containing maps corresponding to the clustering levels for which maps have been requested using parameters k1 and k2. By default (title=NULL), the pdf file is called "Constrained clustering maps.pdf" if the clustering was produced by \code{\link{constrained.clust}}, or "Unconstrained clustering maps.pdf" if it was produced by \code{\link{hclust}}.

}

\author{ Pierre Legendre, Universite de Montreal, 
based on Fortran function hclust.f by F. Murtagh, modified for R by Ross Ihaka and Fritz Leisch.}

\seealso{ \code{\link{hclust}} }

\references{

  Legendre, P. and L. Legendre. 2012. Numerical ecology, 3rd English
  edition. Elsevier Science BV, Amsterdam.

  Murtagh, F. and P. Legendre. Ward’s hierarchical agglomerative clustering method: which algorithms implement Ward’s criterion? (Submitted manuscript).

  Ward, J. H. 1963. Hierarchical grouping to optimize an objective function. 
  J. Amer. Statist. Assoc. 58: 236-244.
}

\examples{ 
###
### Artificial map data from Legendre & Legendre (2012, Fig. 13.26): n = 16

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

# 4. Plot maps of constrained clustering steps from 2 to 8 groups
map.hclust(res.dat, coord.dat)

# Additional operations

# 5. Plot dendrogram (if needed)
# Reversals are caused by the spatial contiguity constraint
plot(res.dat, hang=-1)

# 6. Obtain list of neighbors and plot a map of the triangulation edges
neighbors <- listw2sn(listW)[,1:2]

plot(coord.dat, type='n',asp=1)
title("Delaunay triangulation")
text(coord.dat, labels=as.character(as.matrix(dat)), pos=3)
for(i in 1:nrow(neighbors)) lines(rbind(coord.dat[neighbors[i,1],], coord.dat[neighbors[i,2],]))

# Alternative function to compute matrix links.mat.dat from a list of links
links.mat.dat <- contiguity.mat(neighbors, 16)

# 7. Unconstrained hierarchical clustering and maps, same data
res2.dat <- hclust(D.dat, method="ward")
plot(res2.dat, hang=-1)
map.hclust(res2.dat, coord.dat, 2, 8)

# 8. Another way (but slower) to obtain unconstrained hierarchical clustering:
# Construct a links.mat matrix with connexions everywhere 
# and use it in the constrained clustering function
n <- 16    # Number of objects
links.mat.1 <- matrix(1, n, n)
diag(links.mat.1) <- 0
res.dat.1 <- constrained.clust(D.dat, links.mat.1)

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

# 4. Plot maps of constrained clustering steps from 2 to 10 groups
map.hclust(res.mite, mite.xy, 2, 10, pos=4)

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
cross.fish$AIC

# 4. Plot maps of constrained clustering steps from 2 to 8 groups
coord.fish <- 1:30
map.hclust(res.fish, coord.fish, 2, 8, cex=1.5, pos=3)

}

\keyword{ multivariate }