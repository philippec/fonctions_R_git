\name{contiguity.mat}
\alias{contiguity.mat}
\title{ Constructs a contiguity matrix }

\description{
  Function \code{contiguity.mat} constructs a contiguity matrix from the data read in a file of connexions between points.
}

\usage{
contiguity.mat(data, n)
}

\arguments{
  \item{data}{ A matrix or data frame containing a list of connected pairs of points. }
  
  \item{n}{ Number of objects involved in the \code{data} list. }

} 

\details{ Using the list of connected points in the \code{data} file, the function constructs a square matrix containing 1 for pairs of connected sites and 0 elsewhere.
}

\value{ 
  \item{mat }{A n by n matrix containing 1 for pairs of connected sites and 0 elsewhere. The diagonal values are 0. }
  
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
listW <- nb2listw(tri2nb(coord.dat))

# 2. Obtain the list of neighbors
# Each edge is present twice in the file, e.g. (1 2) and (2 1).
neighbors <- listw2sn(listW)[,1:2]

# Edges can be removed from the list by editing the file
# To eliminate an edge, remove both copies from the file

# 3. Plot a map of the triangulation edges
plot(coord.dat, type='n',asp=1)
title("Delaunay triangulation")
text(coord.dat, labels=1:16, pos=3)
for(i in 1:nrow(neighbors)) lines(rbind(coord.dat[neighbors[i,1],], coord.dat[neighbors[i,2],]))

# Data values can be plotted on the map instead of object numbers
text(coord.dat, labels=as.character(as.matrix(dat)), pos=3)

# 4. Compute matrix links.mat from the list of neighbor links
links.mat.dat <- contiguity.mat(neighbors, 16)

###
### Example of spatial transect or time series

# 1. The transect or time series is represented by a line with equispaced points. The list of connections is created using the rank order of the points.
from.to <- matrix(c(1:29, 2:30), 29, 2)
colnames(from.to) <- c("From", "To")

# 2. Construct the links.mat matrix
links.mat.fish <- contiguity.mat(from.to, 30)
}

\keyword{ multivariate }