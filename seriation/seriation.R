seriation <- function(mat, verbose=FALSE)
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
	
	mat.original <- mat.tmp
	
	# The diagonal elements are removed from the calculations
	diag(mat.tmp) <- NA
	if(verbose) {
		print(mat.tmp)
		cat('\n')
		}

	# 2. Assign weights from n to 1 to the rows
	mult <- c(nb.row:1)
	
	sdev.tmp <- 0.1
	sdevmax <- 0
	sdev.data1 <- 0
	sdev.data2 <- 0
	count <- 0
	arret <- FALSE
	while(sdev.tmp != sdev.data1 & sdev.tmp != sdev.data2 & !arret & count < 1000) {
	# while(count < 6) {
		count <- count+1
		sdev.data2 <- sdev.data1
		sdev.data1 <- sdev.tmp
		
		# 1. Compute column sums, excluding diagonal elements
		sum.col.tmp <- apply(mat.tmp, 2, sum, na.rm=TRUE)
		if(verbose) cat("sum.col.tmp",sum.col.tmp,'\n')
		
		# 3. Multiply table elements by row weights 'mult'
		#    Compute new column sums, excluding diagonal elements
		sum.pr.tmp <- apply((mat.tmp*mult), 2, sum, na.rm=TRUE)
		if(verbose) cat("sum.pr.tmp",sum.pr.tmp,'\n')
		
		# 4. Compute mean weights
		pr.ave.tmp <- sum.pr.tmp/sum.col.tmp
		# Compute standard deviation of mean weights, to be maximized
		sdev.tmp <- sd(pr.ave.tmp)
		if(verbose) cat("sdev.tmp",sdev.tmp,'\n','\n')
		
		# 5. Compute ranks of mean weights (rank 1 to the largest value)
		ranking.tmp <- sort(pr.ave.tmp, decreasing=TRUE, index.return=TRUE)$ix
		# if(verbose) cat("ranking.tmp",ranking.tmp,'\n')
		
		# 6. Reorder the rows/columns in decreasing order of pr.ave for next iteration
		if(sdev.tmp >= sdevmax) {
			sdevmax <- sdev.tmp
			sum.col <- sum.col.tmp
			sum.pr <- sum.pr.tmp
			pr.ave <- pr.ave.tmp
			mat <- mat.tmp        # Save previous matrix in 'mat' before rearrangement
			mrk <- mrk.tmp
			mat.tmp <- mat.tmp[ranking.tmp,ranking.tmp]
			mrk.tmp <- mrk.tmp[ranking.tmp]
			} else { arret <- TRUE }

		if(verbose & !arret) {
			cat('Iteration',(count-1),'sdev =',sdev.tmp,'\n')
			cat('Iteration',count,'ranking.tmp =',ranking.tmp,'\n')
			# cat('mrk.tmp     =',mrk.tmp,'\n')
			print(mat.tmp)
			cat('\n')
			}
	}
	if(count == 1000) stop("No convergence after 1000 iterations")
	diag(mat) <- 1
	return(list(starting.S=mat.original, reordered.S=mat,ranking=mrk, 			
		dat=rbind(sum.col,sum.pr,pr.ave),sdev=sdevmax))
}
