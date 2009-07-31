\name{plotKM}
\alias{plotKM}
\alias{orderingKM}
\alias{pregraphKM}

\title{ Plot the cascadeKM results }
\description{
 This function plots the partitions produced by the \code{\link{cascadeKM}} K-means wrapper function over the set of objects and rearranges the colours for best representation.
}
\usage{
plotKM(wrapres, min.g, max.g, grpmts.plot=TRUE, sortg=FALSE, gridcol=NA)
}
\arguments{
\item{ wrapres }{ This is the output of the \code{\link{cascadeKM}} function. }
\item{ min.g }{ The minimum number of groups to be displayed. }
\item{ max.g }{ The maximum number of groups to be displayed. }
\item{ grpmts.plot }{ Show the plot (TRUE or FALSE). }
\item{ sortg }{ Sort the objects as a function of their group membership to produce a more easily interpretatable graph. See \eqn{Details}. The original object names are kept; they are used as labels in the output table \code{x}, although not in the graph.  If there were no rownames, sequential row numbers are used to keep track of the original order of the objects.  }
\item{ gridcol }{ The colour of the grid lines in the plots. NA, which is the default value, removes the grid lines.}
}

\details{
  Two plots are produced. The graph on the left has the objects in abscissa and the number of groups in ordinate. The groups are represented by colours. The graph on the right shows the values of the criterion (\code{"calinski"} or \code{"ssi"}) for determining the best partition. In a simulation study, Milligan and Cooper (1985) found that the Calinski-Harabasz criterion recovered the correct number of groups the most often. We recommend that criterion because, if the groups are of equal sizes, the maximum value of "calinski" usually indicates the correct number of groups. Another available index is the simple structure index "ssi". Users should not take the indications of these indices literally when the groups are not equal in size and explore the groups corresponding to other values of \eqn{K}. \cr
  The highest value of the criterion is marked in red. Points marked in orange, if any, indicate partitions producing an increase in the criterion value as the number of groups increases; they may represent other interesting partitions. \cr
  If sortg=TRUE, the objects are reordered by the following procedure: (1) a simple matching distance matrix is computed among the objects, based on the table of K-means assignments to groups, from \eqn{K} = \code{min.g} to \eqn{K} = \code{max.g}. (2) A principal coordinate analysis (PCoA, Gower 1966) is computed on the centred distance matrix. (3) The first principal coordinate is used as the new order of the objects in the graph. A simplified algorithm is used to compute the first principal coordinate only, using the iterative algorithm described in Legendre & Legendre (1998, Table 9.10). The full distance matrix among objects is never computed; this evacuates the problem of storing it when the number of objects is large. Distance values are computed as they are needed by the algorithm. 
}
\value{
\item{ x=x }{ Table with the partitions found for different numbers of groups \eqn{K}, from \eqn{K} = \code{min.g} to \eqn{K} = \code{max.g}. }
\item{ best.grps}{ The number of objects in each group (rows) of each partition (columns) from \eqn{K} = \code{min.g} to \eqn{K} = \code{max.g}. }
}
\references{ 
Calinski, T. and J. Harabasz. 1974. A dendrite method for cluster
analysis. Commun. Stat. 3: 1-27. \cr
Gower, J. C. 1966. Some distance properties of latent root and vector methods used in multivariate analysis. Biometrika 53: 325-338. \cr
Legendre, P. & L. Legendre. 1998. Numerical ecology, 2nd English edition.
Elsevier Science BV, Amsterdam. \cr
Milligan, G. W. & M. C. Cooper. 1985. An examination of procedures for determining the number of clusters in a data set. Psychometrika 50: 159-179.
 }
\author{ Sebastien Durand <Sebastien.Durand@UMontreal.ca> and Pierre Legendre <Pierre.Legendre@UMontreal.ca>
}
\seealso{ \code{\link{cascadeKM}}, \code{\link{cIndexKM}}}
\examples{
 # Partitioning a (10 x 10) data matrix of random numbers
 mat <- matrix(runif(100),10,10)
 res <- cascadeKM(mat, 2, 5, iter = 25, criterion = 'calinski') 
 toto <- plotKM(res)
 
 # Partitioning an autocorrelated time series
 vec<-sort(matrix(runif(30),30,1))
 res <- cascadeKM(vec, 2, 5, iter = 25, criterion = 'calinski')
 plotKM(res)
 
 # Add grid lines
 toto <- plotKM(res, gridcol="black")
 
 # Partitioning a larger number of objects
 res <- cascadeKM(matrix(runif(300),300,1), 2, 7,iter = 50, criterion = 'calinski') 
 toto <- plotKM(res)
 
 # Sort the objects based on their group membership to ease interpretation
 toto <- plotKM(res, gridcol=NA, sortg=TRUE)

 # Partitioning a large autocorrelated time series
 # Remove the grid lines
 vec <- sort(matrix(runif(1000),1000,1))
 res <- cascadeKM(vec, 2, 7, iter = 10, criterion = 'calinski')
 toto <- plotKM(res)

 # Use the "ssi" criterion to determine the best partition
 res <- cascadeKM(vec, 2, 12, iter = 20, criterion = "ssi")
 toto <- plotKM(res,4,10)
    
}
\keyword{ cluster }
\keyword{ dplot }
