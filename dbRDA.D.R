dbRDA.D <- function(D, X, nperm=999, compute.eig=FALSE, coord=FALSE)
#
# Compute the dbRDA F-test. The response is represented by a Euclidean or 
# non-Euclidean dissimilarity matrix and X is a matrix of explanatory variables, 
# as in regular RDA. 

# The function uses the computation method of McArdle and Anderson (2001). 
# The F-statistic is obtained without computation of the eigenvalues and 
# eigenvectors, hence no correction has to be made for the negative eigenvalues.
#
# The eigenvalues and eigenvectors are also computed if compute.eig=TRUE.
# The eigenvectors are not scaled to sqrt(eigenvalues), hence they are
# not principal coordinates in the PCoA sense. 
# Optionally, is coord=TRUE, the principal coordinates corresponding to the 
# positive eigenvalues of D are computed.
#
# Arguments --
#
# D : Distance matrix representing the response data. D may be non-Euclidean.
# X : Matrix of explanatory variables for the RDA
# nperm : Number of permutations for the test of significance.
# compute.eig=TRUE if the eigenvalues and eigenvectors should be computed.
# coord=TRUE to compute the principal coordinates corresponding to the
#    positive eigenvalues of D. Requires that compute.eig=TRUE.
#
# Value --
#
# F : F-statistic.
# P.perm : Permutational p-value.
# SS.total : Trace of matrix G, equal to the total sum of squares and to the 
#            sum of the eigenvalues of D.
# values : Eigenvalues (if they are computed, i.e. if compute.eig=TRUE).
# vectors : Eigenvectors (if they are computed, i.e. if compute.eig=TRUE).
# coord : Principal coordinates corresponding to the positive eigenvalues of D.
#
# References --
#
# Legendre, P. and L. Legendre. 2012. Numerical ecology, 3rd English edition. 
# Elsevier Science BV, Amsterdam. 
#
# McArdle, B. H. and M. J. Anderson. 2001. Fitting multivariate models to 
# community data: a comment on distance-based redundancy analysis. 
# Ecology 82: 290-297. 
#
# Example -- Using some sites from the mite data
#
# library(vegan)
# library(PCNM)
# data(mite)
# data(mite.xy)
# mite.BC = vegdist(mite.[10:20,], "bray")   # One negative eigenvalue
# mite.mem = PCNM(dist(mite.xy[10:20,]), dbMEM=TRUE)   # Three MEM are produced
# Load function "dbrda.D"
# res = dbRDA.D(mite.BC,mite.mem$vectors,nperm=999,compute.eig=TRUE,coord=TRUE)
#
# License: GPL-2
# Author:: Pierre Legendre, December 2012
{
	D <- as.matrix(D)
	n <- nrow(D)
	m <- qr(X, tol=1e-6)$rank
	X <- scale(X, center=TRUE, scale=FALSE)   # Centre matrix X
	epsilon <- sqrt(.Machine$double.eps) 
#
# Gower centring, matrix formula. Legendre & Legendre (2012), equation 9.42
	One <- matrix(1,n,n)
	mat <- diag(n) - One/n
	G <- -0.5 * mat %*% (D^2) %*% mat
	trace <- sum(diag(G))
	# LCBD <- diag(G)
#
# Eigenvalue decomposition
	if(compute.eig) {
		eig <- eigen(G, symmetric=TRUE)
		values <- eig$values     # All eigenvalues
		vectors <- eig$vectors   # All eigenvectors, scaled to lengths 1
		if(coord) {
			select <- which(values > epsilon)
			princ.coord <- vectors[,select] %*% diag(sqrt(values[select]))
			} else { princ.coord <- NA }
		} else {
		values <- vectors <- princ.coord <- NA
		}
# Compute the F statistic: McArdle & Anderson (2001), equation 4
# H is the projector matrix, also called the "hat" matrix in the stat literature
	H <- X %*% solve(t(X) %*% X) %*% t(X)
	I.minus.H <- diag(n) - H
	num <- sum(diag(H %*% G %*% H))
	den <- sum(diag(I.minus.H %*% G %*% I.minus.H))
	F <- num/den
# Permutation test of F
nGE=1
for(i in 1:nperm)
	{
	order <- sample(n)
	G.perm <- G[order, order]
	num <- sum(diag(H %*% G.perm %*% H))
	den <- sum(diag(I.minus.H %*% G.perm %*% I.minus.H))
	F.perm <- num/den
	if(F.perm >= F) nGE=nGE+1
	}
P.perm <- nGE/(nperm+1)
#
list(F=F*(n-m-1)/m, P.perm=P.perm, SS.total=trace, values=values, vectors=vectors, coord=princ.coord)
# list(F=F*(n-m-1)/m, P.perm=P.perm, G=G, H=H, SS.total=trace, values=values, vectors=vectors, coord=princ.coord)
}
