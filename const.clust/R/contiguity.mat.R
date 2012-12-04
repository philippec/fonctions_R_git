'contiguity.mat' <- 
	function(data, n)
# This function constructs a contiguity matrix from the data read from
# a file of connexions between points.
#
# Licence: GPL-2
# Author: Pierre Legendre
{
mat <- matrix(0, n, n)
# mat1 <- matrix(NA, n, n)
for(k in 1:nrow(data)) {
	i <- data[k,1]
	if(i > n) stop("Error: 'From' > n")
	j <- data[k,2]
	if(i > n) stop("Error: 'To' > n")
	mat[i,j] <- mat[j,i] <- 1
	# mat1[i,j] <- mat1[j,i] <- 1
	}
mat
}