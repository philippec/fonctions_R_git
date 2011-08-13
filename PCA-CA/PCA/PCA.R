`PCA` <- 
   function(Y, stand=FALSE, cumfit.obj=TRUE, cumfit.var=TRUE)
# 
# Principal component analysis (PCA) with option for variable standardization
#
# stand = FALSE : center by columns only, do not divide by s.d.
# stand = TRUE  : center and standardize (divide by s.d.) by columns
#
# cumfit.obj = TRUE: compute cumulative fit of objects
# cumfit.var = TRUE: compute cumulative fit of variables
#
#          Pierre Legendre, May 2006
{
# Begin internal functions
sq.length <- function(vec) sum(vec^2)
#
'cumul.fit.var' <- function(Y,U2,n,p,k,var.names)
	# Compute the table of "Cumulative fit per variable" 
{
	sp.var <- diag(var(Y))
	res <- matrix(NA,p,k)
	for(i in 1:p) {
	   res[i,] <- cumsum(U2[i,]^2)/sp.var[i]
	   }
	rownames(res) <- var.names
	colnames(res) <- paste("Cum.axis",1:k,sep=".")
	res
}
#
'cumul.fit.obj' <- function(Y,F,n,p,k,obj.names)
	# Compute the table of "Cumulative fit of the objects" 
{
	sq.length <- function(vec) sum(vec^2)
	res <- matrix(NA,n,k)
	for(i in 1:n) {
	   res[i,] <- cumsum(F[i,]^2)/sq.length(F[i,])
	   }
	rownames(res) <- obj.names
	colnames(res) <- paste("Cum.axis",1:k,sep=".")
	res
}
# End internal functions
#
   Y <- as.matrix(Y)
   obj.names <- rownames(Y)
   var.names <- colnames(Y)
   n <- nrow(Y)
   p <- ncol(Y)
   Y.cent <- apply(Y, 2, scale, center=TRUE, scale=stand)
   Y.cov <- cov(Y.cent)
   Y.eig <- eigen(Y.cov)
   k <- length(which(Y.eig$values > 1e-10))
   U  <- Y.eig$vectors[,1:k]
   F  <- Y.cent %*% U
   U2 <- U %*% diag(Y.eig$value[1:k]^(0.5))
   G  <- F %*% diag(Y.eig$value[1:k]^(-0.5))
   rownames(F)  <- obj.names
   rownames(U)  <- var.names
   rownames(G)  <- obj.names
   rownames(U2) <- var.names
   axenames <- paste("Axis",1:k,sep=" ")
   colnames(F)  <- axenames
   colnames(U)  <- axenames
   colnames(G)  <- axenames
   colnames(U2) <- axenames
#
# Fractions of variance
   varY <- sum(diag(Y.cov))
   eigval <- Y.eig$values[1:k]
   relative <- eigval/varY
   rel.cum <- cumsum(relative)
#
   if(cumfit.var) {
      cfit.var <- cumul.fit.var(Y.cent,U2,n,p,k,var.names)
      } else {
      cfit.var <- NULL
      }
#
   if(cumfit.obj) {
      cfit.obj <- cumul.fit.obj(Y.cent,F,n,p,k,obj.names)
      } else {
      cfit.obj <- NULL
      }
#
out <- list(total.var=varY, eigenvalues=eigval, rel.eigen=relative, 
       rel.cum.eigen=rel.cum, U=U, F=F, U2=U2, G=G, 
       cumulative.fit.var=cfit.var, cumulative.fit.obj=cfit.obj, stand=stand, 
       obj.names=obj.names, var.names=var.names, call=match.call() )
class(out) <- "PCA"
out
}

`print.PCA` <-
    function(x, kk=5, ...)
{
    cat("\nPrincipal Component Analysis\n")
    cat("\nCall:\n")
    cat(deparse(x$call),'\n')
    if(x$stand) cat("\nThe data have been centred and standardized by column",'\n')
    cat("\nTotal variance in matrix Y: ",x$total.var,'\n')
    cat("\nEigenvalues",'\n')
    cat(x$eigenvalues,'\n')
    cat("\nRelative eigenvalues",'\n')
    cat(x$rel.eigen,'\n')
    cat("\nCumulative relative eigenvalues",'\n')
    cat(x$rel.cum.eigen,'\n')
    kk <- min(length(x$eigenvalues), kk)
    if(!is.null(x$cumulative.fit.var)) {
       cat("\nCumulative fit per variable (",kk,"axes)",'\n')
       print(x$cumulative.fit.var[,1:kk])
       }
    if(!is.null(x$cumulative.fit.obj)) {
       cat("\nCumulative fit of the objects (",kk,"axes)",'\n')
       print(x$cumulative.fit.obj[,1:kk])
       }
    cat('\n')
    invisible(x) 
}

`biplot.PCA` <-
    function(x, scaling=1, plot.axes=c(1,2), color.obj="black", color.var="red", ...)
# scaling = 1 : preserves Euclidean distances among the objects
# scaling = 2 : preserves correlations among the variables
{
    #### Internal function
	larger.frame <- function(mat, percent=0.07)
	# Produce an object plot 10% larger than strictly necessary
	{
	range.mat = apply(mat,2,range)
	z <- apply(range.mat, 2, function(x) x[2]-x[1])
	range.mat[1,]=range.mat[1,]-z*percent
	range.mat[2,]=range.mat[2,]+z*percent
	range.mat
	}
	####
	
	if(length(x$eigenvalues) < 2) stop("There is a single eigenvalue. No plot can be produced.")
	if(length(which(scaling == c(1,2))) == 0) stop("Scaling must be 1 or 2")

	par(mai = c(1.0, 0.75, 1.0, 0.5))

	if(scaling == 1) {

	# Distance biplot, scaling type = 1: plot F for objects, U for variables
	# This projection preserves the Euclidean distances among the objects
	
	lf.F = larger.frame(x$F[,plot.axes])
	biplot(x$F[,plot.axes],x$U[,plot.axes],col=c(color.obj,color.var), xlim=lf.F[,1], ylim=lf.F[,2], arrow.len=0.05, asp=1)
	title(main = c("PCA biplot","scaling type 1"), family="serif", line=3)

	} else {

	# Correlation biplot, scaling type = 2: plot G for objects, U2 for variables
	# This projection preserves the correlation among the variables
   
	lf.G = larger.frame(x$G[,plot.axes])
	biplot(x$G[,plot.axes],x$U2[,plot.axes],col=c(color.obj,color.var), xlim=lf.G[,1], ylim=lf.G[,2], arrow.len=0.05, asp=1)
	title(main = c("PCA biplot","scaling type 2"), family="serif", line=3)

	}
invisible()
}

# mite.hel = decostand(mite, "hel")
# mite.hel.D = dist(mite.hel)
# mite.correlog = mantel.correlog(mite.hel.D, XY=mite.xy, nperm=99, cutoff=FALSE)
# mite.correlog
# mite.correlog$mantel.res
# plot(mite.correlog)
