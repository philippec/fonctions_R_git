'constrained.clust' <- 
	function(D, links.mat=NULL, method="ward.D2", beta=-0.25, verbose=FALSE, target=1)
# R-coded version of the hclust.f function in Fortran
# modified to produce spatially-constrained clustering solutions as per
# Legendre, P. and L. Legendre, 1998. Numerical ecology. Elsevier (p. 758).
#
# Parameters -
# D : dissimilarity matrix
# links.mat : contiguity matrix among the points (spatial connexions, 0-1)
# methods : "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), 
#           "mcquitty" (= WPGMA), "centroid" (= UPGMC), "median" (= WPGMC)",
#           "flexible"
#
# Parameters in the function
# n   : number of objects
# len = n*(n-1)/2
# diss(len) : dissimilarity vector, evolving as clustering proceeds
# method {1:7} : clustering method
# membr(n) : vector of cluster cardinalities, length = n
# nn(n) : vector of current nearest neighbour of each object, length = n
# crit : D-value associated with nearest neighbour
# flag(n) : vector of boolean indicator of agglomerable objects or clusters
#
# License: GPL-2 
# Author:: Pierre Legendre
# based on Fortran function hclust.f by F. Murtagh, modified for R by R. Hiaka

{
### Internal functions: nearest neighbour (NN) chain agglomeration
#
# Find positions of values in upper-triangular D turned into a vector. 
# Always use this function with i < j.
ioffset <- function(n,i,j)  { j + (i-1)*n - (i*(i+1))/2 }
#
LW <- function(method, diss, i, i2, j2, ind1, ind2, D12, membr, alpha, beta)
# Lance & Williams agglomerative algorithm as per Fionn Murtagh
{
if(method <= 2) {# ward.D and ward.D2
	temp = (membr[i2]+membr[i])*diss[ind1] +
	(membr[j2]+membr[i])*diss[ind2] - membr[i]*D12
	res = temp / (membr[i2]+membr[j2]+membr[i])

	} else if(method==3) { # Single linkage
	res = min(diss[ind1], diss[ind2])
	
	} else if(method==4) { # Complete linkage
	res = max(diss[ind1], diss[ind2])
	
	} else if(method==5) { # UPGMA ("average")
	res = (membr[i2]*diss[ind1] + membr[j2]*diss[ind2])/
	(membr[i2]+membr[j2])
     
	} else if(method==6) { # WPGMA ("mcquitty")
	res = 0.5*diss[ind1] + 0.5*diss[ind2]
	
	} else if(method==7) { # UPGMC ("centroid")
	temp = (membr[i2]*diss[ind1] + membr[j2]*diss[ind2] -
	membr[i2]*membr[j2]*D12/(membr[i2] + membr[j2]))
	res = temp / (membr[i2]+membr[j2])
	
	} else if(method==8) { # WPGMC ("median")
	res = 0.5*diss[ind1] + 0.5*diss[ind2] - 0.25*D12
	
	} else if(method==9) { # Flexible
	res = alpha*diss[ind1] + alpha*diss[ind2] + beta*D12
	
	}
res
}
### End internal functions

	# Preliminaries
    METHODS <- c("ward.D", "ward.D2", "single", "complete", "average",  
        "mcquitty", "centroid", "median", "flexible")
    method <- pmatch(method, METHODS)
    if (is.na(method)) stop("invalid clustering method")
    n <- as.integer(attr(D, "Size"))
    if (is.null(n)) stop("invalid dissimilarities")
    if (n < 2) stop("must have n >= 2 objects to cluster")
    len <- as.integer(n * (n - 1)/2)
    if (length(D) != len) 
    	(if (length(D) < len) stop
		else warning) ("D matrix of improper length")
	if(method==2) D = D^2
	if((beta < -1) | (beta>= 1)) stop("Flexible clustering: (-1<= beta < 1)")
    alpha = (1-beta)/2

    # If there is no links matrix ...
    if(is.null(links.mat)) {
    	links.mat <- matrix(1,n,n)
    	diag(links.mat) <- 0
    	}
    if((nrow(links.mat) != n) | (ncol(links.mat) != n)) 
    	stop("D and links.mat do not match in numbers of objects")

	# Initialization    
    membr <- rep(1, n)
    flag <- as.logical(rep(1, n))
    ncl <- n
    disnn = vector(mode="numeric", length=n)
    nn    = vector(mode="integer", length=n)
    ia    = vector(mode="integer", length=n)
    ib    = vector(mode="integer", length=n)
    crit  = vector(mode="numeric", length=n)
    
    # Create vector from upper-triangular D by rows, 
    # or lower-triangular by columns (identical vector)
    diss0  = as.vector(D)
    
    # Hadamard product of D by links.mat
    # Use strung-out matrices, i.e. transformed into vectors
    #
    # Matrix "mat" has NA for unconnected objects, instead of 0 in links.mat.
    # This will allow the algorithm to recognize and handle D values of 0:
    # connected objects that have D = 0 will produce a value 0 in "diss",
    # whereas unconnected objects will produce a value NA.
    mat <- matrix(NA,n,n)
	mat[links.mat==1] <- 1
	# Hadamard product of strung-out matrices
    diss = as.vector(as.dist(mat)) * diss0  
    
    # Carry out an agglomeration - First create list of nearest neighbours.
	# Note: nn contains the nearest neighbour of each object TO THE RIGHT of i.
	# disnn contains the distance of each object to its nearest neighbour.
	
	# Create the list of nearest neighbours
	for(i in 1:(n-1)) {
		Dmin <- .Machine$double.xmax
		for(j in (i+1):n) {
			ind <- ioffset(n,i,j)
			if(!is.na(diss[ind]) & (diss[ind] < Dmin)) {
				Dmin <- diss[ind]
				jm <- j
				}
			}
		nn[i] <- jm
		disnn[i] <- Dmin
		}
    	
###
### Begin clustering loop
	while(ncl > 1) {
		if(verbose) {
			cat('Number of groups (ncl) =',ncl,'\n')
			cat('nn[target] =',nn[target],'\n')
			cat('disnn[target] =',disnn[target],'\n')
			vec1 = ioffset(n,(1:(target-1)),target)
			vec2 = ioffset(n,target,(target+1):n)
			cat('diss0[target] =',diss0[vec1],0,diss0[vec2],'\n')
			cat('links.mat[target,] =',links.mat[target,],'\n')
			cat('flag =',flag,'\n')
			cat('\n')
			}

    	# Determine next pair to cluster using list of nearest neighbours
		Dmin <- .Machine$double.xmax
		for(i in 1:(n-1)) {
			if(flag[i] & (!is.na(disnn[i])) & (disnn[i] < Dmin)) {
				Dmin <- disnn[i]
				im = i
				jm = nn[i]
				}
			}

       	ncl <- ncl-1

		# Agglomerate
		i2 = min(im,jm)    # First object in the clustering pair
		j2 = max(im,jm)    # Second object in the clustering pair
		ia[n-ncl] = i2     # Vector of first objects in pairs
		ib[n-ncl] = j2     # Vector of second objects in pairs
		crit[n-ncl] = Dmin
		flag[j2] = FALSE
		
		# Update diss0 matrix for new cluster
		for(i in 1:n) {
			if(flag[i] & (i!=i2)) {
				if(i2<i) {ind1=ioffset(n,i2,i)} else {ind1=ioffset(n,i,i2)}
				if(j2<i) {ind2=ioffset(n,j2,i)} else {ind2=ioffset(n,i,j2)}
				ind3 = ioffset(n,i2,j2)
				D12 = diss0[ind3]
				#
				diss0[ind1] <- 
				LW(method, diss0,i,i2,j2,ind1,ind2,D12,membr,alpha,beta)
				}
			}
		membr[i2] <- membr[i2]+membr[j2]
		
		# Update links.mat: fuse rows and columns of i2 and j2
		vec <- links.mat[i2,] + links.mat[j2,]
		vec[vec > 1] = 1
		links.mat[i2,] <- as.numeric(vec)
		links.mat[,i2] <- as.numeric(vec)
		mat <- matrix(NA,n,n)
		mat[links.mat > 0] <- 1
		
		# Hadamard product of strung-out matrices
    	diss = as.vector(as.dist(mat)) * diss0  
		
		# Update list of nearest neighbours
		for(i in 1:(n-1)) {
			if(flag[i]) {
				Dmin <- .Machine$double.xmax
				for(j in (i+1):n) {
					if(flag[j]) {
						ind = ioffset(n,i,j)
						if(!is.na(diss[ind]) & (diss[ind] < Dmin)) {
							Dmin <- diss[ind]
							jj <- j
							}
						}
					}
				nn[i] <- jj
				disnn[i] <- Dmin
				}
			}
				
    	} # End While
### End clustering loop
###
		if(verbose) {
			target = 3
			cat('ncl =',n:2,'\n')
			cat('First object in pairs (ia) =',ia,'\n')
			cat('Second object in pairs (ib) =',ib,'\n')
			cat('crit =',crit,'\n')
			cat('\n')
			}
	if(method==2) crit <- sqrt(crit)

### Prepare results for output
    hcl <- list(n=n, len=len, method=method, ia=ia, ib=ib, crit=crit, 
    		members=membr, nn=nn, disnn=disnn, flag=flag, diss=diss)

    hcass <- .Fortran("hcass2", n = as.integer(n), ia = as.integer(hcl$ia), 
        ib = as.integer(hcl$ib), order = integer(n), iia = integer(n), 
        iib = integer(n))

    tree <- list(merge = cbind(hcass$iia[1L:(n - 1)], hcass$iib[1L:(n - 
        1)]), height = hcl$crit[1L:(n - 1)], order = hcass$order, 
        labels = attr(D, "Labels"), method = METHODS[method], 
        call = match.call(), dist.method = attr(D, "method"), hcl=hcl)
    class(tree) <- "hclust"
    tree
}
