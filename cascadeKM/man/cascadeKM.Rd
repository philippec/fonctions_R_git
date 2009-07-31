\name{cascadeKM}
\alias{cascadeKM}

\title{K-means partitioning using a range of values of K}
\description{
 This function is a wrapper for the \code{kmeans} function of the \code{stats} package. It creates several partitions forming a cascade from a small to a large number of groups.
}

\usage{
 cascadeKM(data, nb_inf_gr, nb_sup_gr, iter = 100, criterion = "calinski")
}

\arguments{
	\item{ data }{ The data matrix. The objects are the rows.}
	\item{ nb_inf_gr }{ The number of groups for the partition with the 
	smallest number of groups of the cascade (min).}
	\item{ nb_sup_gr }{ The number of groups for the partition with the largest 	
	number of groups of the cascade (max).}
	\item{ iter }{ The number of random starting configurations for each value
	of \eqn{K}.}
	\item{ criterion }{ The criterion that will be used to select the best
	partition. The default value is "calinski", which refers to the
	Calinski-Harabasz (1974) criterion. The simple structure index, "ssi", is
	also available. See the \code{\link{cIndexKM}} help file for interpretation
	of these indices. }
}
\details{
  This is a beta version.
}
\value{
	\item{ partition }{ Table with the partitions found for different numbers 
	of groups \eqn{K}, from \eqn{K} = \code{nb_inf_gr} to \eqn{K} = \code{nb_sup_gr}. }
	\item{ results }{ Values of the criterion to select the best partition. }
	\item{ criterion }{ The name of the criterion used. }
	\item{ size }{ The number of objects found in each group, for all 
	partitions (columns). }
}
\references{ Calinski, T. and J. Harabasz. 1974. A dendrite method for cluster
analysis. Commun. Stat. 3: 1-27. }

\author{ Marie-Helene Ouellette <Marie-Helene.Ouellette@UMontreal.ca>,
Sebastien Durand <Sebastien.Durand@UMontreal.ca> and
Pierre Legendre <Pierre.Legendre@UMontreal.ca>
}

\seealso{\code{\link{kmeans}}, \code{\link{plotKM}}, \code{\link{cIndexKM}}}
}
\examples{
 # Partitioning a (10 x 10) data matrix of random numbers
 mat <- matrix(runif(100),10,10)
 res <- cascadeKM(mat, 2, 5, iter = 25, criterion = 'calinski') 
 toto <- plotKM(res)
 
 # Partitioning an autocorrelated time series
 vec <- sort(matrix(runif(30),30,1))
 res <- cascadeKM(vec, 2, 5, iter = 25, criterion = 'calinski')
 toto <- plotKM(res)
 
 # Partitioning a large autocorrelated time series
 # Note that we remove the grid lines
 vec <- sort(matrix(runif(1000),1000,1))
 res <- cascadeKM(vec, 2, 7, iter = 10, criterion = 'calinski')
 toto <- plotKM(res, gridcol=NA)
 
}
\keyword{ cluster }
