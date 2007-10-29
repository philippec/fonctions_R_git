pcnm.all <- function(matdist,thresh=give.thresh(matdist))
#
# Compute the PCNM eigenfunctions corresponding to all eigenvalues (+, 0, -).
# Input file: distance matrix produced by the function "dist".
# Computation of the threshold will require a function of the library "ape".
# 
# Original PCNM function: Stephane Dray, November 11, 2004
# The present version: Pierre Legendre, August 19, 2007
{
    a <- system.time({
	cat("Truncation level =",thresh,'\n')
    matdist <- as.matrix(matdist)
    mattrunc <- matdist
    for (i in 1:dim(matdist)[1]){
        for (j in 1:dim(matdist)[2]){
            mattrunc[i,j] <- ifelse(matdist[i,j]>thresh, 4*thresh,matdist[i,j])
        }
    }
    mypcnm.all <- PCoA.all(mattrunc)
    })
    a[3] <- sprintf("%2f",a[3])
    cat("Time to compute PCNMs =",a[3]," sec",'\n')
    return(mypcnm.all)
}

PCoA.all <- function(D)
# Principal coordinate decomposition of a square distance matrix D
# Get the eigenvectors corresponding to all eigenvalues, positive and negative
# Pierre Legendre, 2005, 2007
{
   DD=as.matrix(D)
   n <- nrow(DD)

# Centring by the matrix formula
   One <- matrix(1,n,n)
   mat <- diag(n) - One/n
   Dpr2 <- -0.5 * mat %*% (DD^2) %*% mat
   trace <- sum(diag(Dpr2))

# Eigenvalue decomposition
   toto <- eigen(Dpr2)
   rel.values <- toto$values/trace
         
   list(eigen.values=toto$values, rel.eigen.values=rel.values, eigen.vectors=toto$vectors, trace=trace)
}

"give.thresh" <- function(distxy)
# This function computes a truncation threshold, such that all points remain
# connected along a minimum spanning tree (mst).
	{
	library(ape)
    require(ape)
    a <- system.time({
    spanning=mst(distxy)
    })
    a[3] <- sprintf("%2f",a[3])
    cat("Time to find the threshold =",a[3]," sec",'\n')
    return(max(spanning*as.matrix(distxy)))
    }
