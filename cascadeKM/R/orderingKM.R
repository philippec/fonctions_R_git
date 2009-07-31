"orderingKM"<-function(mat) {

# INPUT :	
# mat		(n x k): n the objects and k the descriptors
# 			This matrix must be integers and numeric
#			And must not be binairy, it is the partition matrix
#			output by cascadeKM
# OUTPUT :  Ordered matrix
# Note:

	#Check up
	if(!is.matrix(mat)) stop("'mat' must be a matrix!")
	if(!is.numeric(mat)) stop("'mat' must be numeric!")
	if(any(is.na(mat))) stop("'NA' value was found in the matrix!")
	if(any(is.infinite(mat))) stop("'Inf' value was found in the matrix!")
	nb.desc=ncol(mat)
	nb.obj=nrow(mat)
	
	scores<-rep(0.0,nb.obj)
#    dyn.load("~/STATSR/_Packages/cascadeKM/src/orderingKM.so")
#	scores<-as.vector(.Fortran("orderdata",as.integer(mat), 
#	as.integer(nb.obj), as.integer(nb.desc), sc=as.double(scores))$sc)
	scores<-as.vector(.Fortran("orderdata",as.integer(mat), 
	as.integer(nb.obj), as.integer(nb.desc), sc=as.double(scores),
	PACKAGE="cascadeKM")$sc)
	scores <-sort(scores, index.return=TRUE)$ix
	mat<-mat[scores,]
	return(mat)
}
