\name{cIndexKM}
\alias{cIndexKM}

\title{Clustering Indices}
\description{
    Object \code{y} is the result of a clustering algorithm such
    as \code{"cclust"}.
    \code{\link{cIndexKM}} calculates the values of several clustering
    indices. The indices can all be independenly used
    to determine the correct number of clusters in a data set.
}
\usage{ cIndexKM ( y, x, index = "all" )
}

\arguments{
\item{y}{Object of class \code{"kmeans"} returned by a clustering algorithm
such as \code{\link{kmeans}}}
\item{x}{Data matrix where columns correspond to variables and rows to
observations}
\item{index}{The available indices are: \code{"calinski"} and \code{"ssi"}. 
Type \code{"all"} to obtain both indices. 
Abbreviations of these names are also accepted.}
}

\details{
	This function is based on the clustIndex function. 
	Some of the criteria were removed from this version because computation 
	errors were generated when only one object was found in a group.\cr
	The default value is "calinski", 
	which refers to the well-known Calinski-Harabasz (1974) criterion. 
	The other available index is the simple structure index "ssi".
	In the case of groups of equal sizes, "calinski" is generally a good
	criterion to indicate the correct number of groups. Users should not
	take its indications literally when the groups are not equal in size.
    Type "all" to obtain both indices.\cr
    \describe{
    	\item \bold{calinski}:
    	\eqn{(SSB/(K-1))/(SSW/(n-K))}, where \eqn{n} is the
    	number of data points and \eqn{K} is the number of clusters.
    	\eqn{SSW} is the sum of squares within the clusters while
    	\eqn{SSB} is the sum of squares among the clusters. This index
    	is simply an \eqn{F} analysis-of-variance statistic.
    	}
    \describe{
    	\item \bold{ssi}:
    	this "Simple Structure Index" 
    	combines three elements which influence the interpretability of a
    	solution, i.e., the maximum difference of each variable between the
    	clusters, the sizes of the most contrasting clusters and the
    	deviation of a variable in the cluster centers compared to its
    	overall mean. These three elements are multiplicatively combined and
    	normalized to give a value between \eqn{0} and \eqn{1}.
    	}
}

\value{
    Returns a vector with the index values. The maximum value of these indices
    is supposed to indicate the best partition. These indices work best with
    groups of equal sizes. When the groups are not of equal sizes, one should 
    not put too much faith in the maximum of these indices, and also explore the
    groups corresponding to other values of \eqn{K}.
}
\references{Andreas Weingessel, Evgenia Dimitriadou and Sara Dolnicar,
    An Examination Of Indexes For Determining The Number
    Of Clusters In Binary Data Sets,\cr
    \url{http://www.wu-wien.ac.at/am/wp99.htm#29}\cr
    and the references therein.
}
\author{Based on Evgenia Dimitriadou and Andreas Weingessel documentation.
Function modified by Sebastien Durand and Pierre Legendre 
from the original "cIndexKM" in "cclust" to accommodate a single response variable as well as singleton groups.
}

\seealso{\code{\link{kmeans}}, \code{\link{cascadeKM}}, \code{\link{plotKM}}}

\examples{
# a 2-dimensional example
x <- rbind(matrix(rnorm(100, sd = 0.3), ncol = 2),
           matrix(rnorm(100, mean = 1, sd = 0.3), ncol = 2))
colnames(x) <- c("x", "y")
(cl <- kmeans(x, 2))
resultindexes <- cIndexKM(cl,x, index="all")
resultindexes   
}
\keyword{ cluster }

