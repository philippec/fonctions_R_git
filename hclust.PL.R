hclust.PL <- function(D, method="ward.D2", beta=-0.25)
#
# R-coded version of the Fortran function hclust.f, originally written by Fionn
# Murtagh (1986) and modified for R by Ross Ihaka (1996), Fritz Leisch (2000) 
# and Martin Maechler (2001).
#
# Parameters -
# D : dissimilarity matrix
# methods : "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), 
#           "mcquitty" (= WPGMA), "centroid" (= UPGMC), "median" (= WPGMC)",
#           "flexible"
#     Notes: "ward.D" is the method called "ward" in function hclust of stats.
#     "ward.D2" is the modified form that preserves Ward's (1963) criterion.
#     "flexible" is beta-flexible clustering (Lance & Williams l966, 1967).
#     For use, refer to the documentation file of hclust {stats}.
#
# Parameters used in the original Fortran subroutine hclust.f --
# n   : number of objects
# len : n*(n-1)/2
# DISS(len) : dissimilarity vector, evolving as clustering proceeds
# IOPT (called method in the main function) : clustering method
# MEMBR[n) : cluster cardinalities
# NN(n) : current nearest neighbour
# CRIT : D-value associated with nearest neighbour
# FLAG(n) : boolean indicator of agglomerable objects or clusters
#
# Value -
# out$hcl : output of this function
# out$hcl.stats : output of hclust() of {stats}  # Dis-implemented
#
# Author: Pierre Legendre
# based on Fortran subroutine hclust.f by Fionn Murtagh (1986), 
# modified for R by Ross Hiaka (1996), Fritz Leisch (2000) 
# and Martin Maechler (2001).

{
### Internal functions -- 
### LW() is Lance & Williams agglomeration
#
LW <- function(IOPT, DISS, k, I2, J2, ind1, ind2, D12, MEMBR, alpha, beta)
# Lance & Williams agglomerative algorithm as per Fionn Murtagh (1986)
{
if(IOPT <= 2) {# ward.D and ward.D2
	temp = (MEMBR[I2]+MEMBR[k])*DISS[ind1] +
	(MEMBR[J2] + MEMBR[k])*DISS[ind2] - MEMBR[k]*D12
	res = temp / (MEMBR[I2]+MEMBR[J2]+MEMBR[k])

	} else if(IOPT==3) { # Single linkage
	res = min(DISS[ind1], DISS[ind2])
	
	} else if(IOPT==4) { # Complete linkage
	res = max(DISS[ind1], DISS[ind2])
	
	} else if(IOPT==5) { # UPGMA ("average")
	res = (MEMBR[I2]*DISS[ind1] + MEMBR[J2]*DISS[ind2])/
	(MEMBR[I2]+MEMBR[J2])
     
	} else if(IOPT==6) { # WPGMA ("mcquitty")
	res = 0.5*DISS[ind1] + 0.5*DISS[ind2]
	
	} else if(IOPT==7) { # WPGMC ("centroid")
	temp = (MEMBR[I2]*DISS[ind1] + MEMBR[J2]*DISS[ind2] -
	MEMBR[I2]*MEMBR[J2]*D12/(MEMBR[I2] + MEMBR[J2]))
	res = temp / (MEMBR[I2]+MEMBR[J2])
	
	} else if(IOPT==8) { # UPGMC ("median")
	res = 0.5*DISS[ind1] + 0.5*DISS[ind2] - 0.25*D12
	
	} else if(IOPT==9) { # Flexible
	res = alpha*DISS[ind1] + alpha*DISS[ind2] + beta*D12
	
	}
res
}
# Find positions of values in upper-triangular D turned into a vector; 
# always use with i < j
ioffset <- function(n,i,j)  { j + (i-1)*n - (i*(i+1))/2 }
#
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

	# Initialization    
    MEMBR <- rep(1, n)
    FLAG <- as.logical(rep(1, n))
    NCL <- n
    DISNN = vector(mode="numeric", length=n)
    NN = vector(mode="integer", length=n)
    IA = vector(mode="integer", length=n)
    IB = vector(mode="integer", length=n)
    CRIT = vector(mode="numeric", length=n)

    # Vector from upper-triangular D by rows, or lower-triangular by columns
    DISS  = as.vector(D)
    
    # Carry out an agglomeration - first create list of NNs
	# Note NN and DISNN are the nearest neighbour and its distance
	# TO THE RIGHT of I.

	for(i in 1:(n-1)) {
		DMIN <- .Machine$double.xmax
		for(j in (i+1):n) {
			ind = ioffset(n,i,j)
			if(DISS[ind] < DMIN) {
				DMIN <- DISS[ind]
				jm <- j
				}
			}
		NN[i] <- jm        # Current nearest neighbour
		DISNN[i] <- DMIN   # Distance to nearest neighbour
		}
    	
###
### Begin clustering loop
	while(NCL > 1) {
    	# Determine next pair to cluster using list of nearest neighbours

		DMIN <- .Machine$double.xmax
		for(i in 1:(n-1)) {
			if(FLAG[i]) {
				if(DISNN[i] < DMIN) {
					DMIN <- DISNN[i]
					IM = i
					JM = NN[i]
					}
				}
			}

       	NCL <- NCL-1

		# Agglomeration
		I2 = min(IM,JM)   # First object in the clustering pair
		J2 = max(IM,JM)   # Second object in the clustering pair
		IA[n-NCL] = I2    # Vector of first objects in pairs
		IB[n-NCL] = J2    # Vector of second objects in pairs
		CRIT[n-NCL] = DMIN
		FLAG[J2] = FALSE
		
		# Update D matrix for new cluster
		for(k in 1:n) {
			if(FLAG[k] & (k!=I2)) {
				if(I2<k) {ind1=ioffset(n,I2,k)} else {ind1=ioffset(n,k,I2)}
				if(J2<k) {ind2=ioffset(n,J2,k)} else {ind2=ioffset(n,k,J2)}
				ind3=ioffset(n,I2,J2)
				D12 = DISS[ind3]
				#
				DISS[ind1] = 
				LW(method, DISS, k, I2, J2, ind1, ind2, D12, MEMBR, alpha, beta)
				}
			}
		MEMBR[I2] <- MEMBR[I2]+MEMBR[J2]
		
		# Update list of nearest neighbours
		for(i in 1:(n-1)) {
			if(FLAG[i]) {
				DMIN <- .Machine$double.xmax
				for(j in (i+1):n) {
					if(FLAG[j]) {
						ind = ioffset(n,i,j)
						if(DISS[ind] < DMIN) {
							DMIN <- DISS[ind]
							jj <- j
							}
						}
					}
				NN[i] <- jj
				DISNN[i] <- DMIN
				}
			}
    	} # End While

	if(method==2) CRIT <- sqrt(CRIT)

    hcl <- list(n=n, len=len, method=method, ia=IA, ib=IB, crit=CRIT, 
    		members=MEMBR, nn=NN, disnn=DISNN, flag=FLAG, diss=DISS)

#	Call function hclust.f of {stats}
#	members <- rep(1, n)
#	hcl.stats <- .Fortran("hclust", n=n, len=len, method=as.integer(method), 
#        ia = integer(n), ib = integer(n), crit = double(n), 
#        members = as.double(members), 
#        nn = integer(n), disnn = double(n), flag = logical(n), 
#        diss = as.double(D), PACKAGE = "stats")

    hcass <- .Fortran("hcass2", n = as.integer(n), ia = as.integer(hcl$ia), 
        ib = as.integer(hcl$ib), order = integer(n), iia = integer(n), 
        iib = integer(n), PACKAGE = "stats")

    tree <- list(merge = cbind(hcass$iia[1L:(n - 1)], hcass$iib[1L:(n - 
        1)]), height = hcl$crit[1L:(n - 1)], order = hcass$order, 
        labels = attr(D, "Labels"), method = METHODS[method], 
        call = match.call(), dist.method = attr(D, "method"), 
        hcl=hcl)
    class(tree) <- "hclust"
    tree
}
