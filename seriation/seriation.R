seriation <- function(mat)
{	
	# Check: similarity or distance matrix
	mat.tmp <- as.matrix(mat)
	if(!is.numeric(mat.tmp)) stop("Data matrix 'mat' must be numeric")
	nb.col <- ncol(mat.tmp)
	nb.row <- nrow(mat.tmp)
	if(nb.col != nb.row) stop("Data matrix 'mat' must be a square matrix")
	if(min(mat.tmp) < 0) stop("Data matrix 'mat' contains negative values")
	
	mrk.tmp <- 1:nb.col
	if(is.null(colnames(mat.tmp))) {
	   colnames(mat.tmp) <- 1:nb.col
	   rownames(mat.tmp) <- 1:nb.col
	   }
	# If a distance matrix is provided, convert it to similarity
	maxmat <- max(mat.tmp)
	if(mat.tmp[1,1] == 0) {
		if(maxmat <= 1) {
		mat.tmp <- 1-mat.tmp
		} else {
		mat.tmp <- 1-(mat.tmp/maxmat)
		}
		}
	
	# The diagonal elements are removed from the calculations
	diag(mat.tmp) <- NA
	
	mult <- c(nb.row:1)
	sdev.tmp <- 0.1
	sdevmax <- 0
	sdev.data1 <- 0
	sdev.data2 <- 0
	count <- 0
	while(sdev.tmp != sdev.data1 & sdev.tmp != sdev.data2 & count < 1000) {
		count <- count+1
		sdev.data2 <- sdev.data1
		sdev.data1 <- sdev.tmp
		# Compute column sums, excluding diagonal elements
		sum.col.tmp <- apply(mat.tmp, 2, sum, na.rm=TRUE)
		# Multiply table elements by row weights 'mult'
		# Compute new column sums, excluding diagonal elements
		sum.pr.tmp <- apply((mat.tmp*mult), 2, sum, na.rm=TRUE)
		# Compute mean weights
		pr.ave.tmp <- sum.pr.tmp/sum.col.tmp
		# Compute standard deviation of mean weights, to be maximized
		sdev.tmp <- sd(pr.ave.tmp)
		# Compute ranks of mean weights (rank 1 to the largest value)
		ranking.tmp <- sort(pr.ave.tmp, decreasing=TRUE, index.return=TRUE)$ix	
		# Reorder the rows/columns in decreasing order of pr.ave for next iteration
		if(sdev.tmp >= sdevmax) {
		   sdevmax <- sdev.tmp
		   sum.col <- sum.col.tmp
		   sum.pr <- sum.pr.tmp
		   pr.ave <- pr.ave.tmp
		   mat <- mat.tmp
		   mrk <- mrk.tmp
		   }
		mat.tmp <- mat.tmp[ranking.tmp,ranking.tmp]
		mrk.tmp <- mrk.tmp[ranking.tmp]
	}
	if(count == 1000) stop("No convergence after 1000 iterations")
	diag(mat) <- 1
	return(list(mat=mat,ranking=mrk, 			
		dat=rbind(sum.col,sum.pr,pr.ave),sdev=sdevmax))
}
