`PCA` <- 
   function(Y, stand=FALSE, scaling=1, color.obj="black", color.var="red")
# 
# Principal component analysis (PCA) with options for scalings and output
#
# stand = FALSE : center by columns only, do not divide by s.d.
# stand = TRUE  : center and standardize (divide by s.d.) by columns
# scaling = 1 : preserves Euclidean distances among the objects
# scaling = 2 : preserves correlations among the variables
#
#          Pierre Legendre, May 2006
{
   Y = as.matrix(Y)
   if(length(which(scaling == c(1,2))) == 0) stop("Scaling must be 1 or 2")
   obj.names = rownames(Y)
   var.names = colnames(Y)
   size = dim(Y)
   Y.cent = apply(Y, 2, scale, center=TRUE, scale=stand)
   Y.cov = cov(Y.cent)
   Y.eig = eigen(Y.cov)
   k = length(which(Y.eig$values > 1e-10))
   U  = Y.eig$vectors[,1:k]
   F  = Y.cent %*% U
   U2 = U %*% diag(Y.eig$value[1:k]^(0.5))
   G  = F %*% diag(Y.eig$value[1:k]^(-0.5))
   rownames(F)  = obj.names
   rownames(U)  = var.names
   rownames(G)  = obj.names
   rownames(U2) = var.names
#
# Fractions of variance
   varY = sum(diag(Y.cov))
   eigval = Y.eig$values[1:k]
   relative = eigval/varY
   rel.cum = vector(length=k)
   rel.cum[1] = relative[1]
   for(kk in 2:k) { rel.cum[kk] = rel.cum[kk-1] + relative[kk] }
#
out <- list(total.var=varY, eigenvalues=eigval, rel.eigen=relative, 
       rel.cum.eigen=rel.cum, U=U, F=F, U2=U2, G=G, stand=stand, 
       scaling=scaling, obj.names=obj.names, var.names=var.names,
       color.obj=color.obj, color.var=color.var, call=match.call() )
class(out) <- "PCA"
out
}

`print.PCA` <-
    function(x, ...)
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
    invisible(x) 
}

`biplot.PCA` <-
    function(x, ...)
{
if(length(x$eigenvalues) < 2) stop("There is a single eigenvalue. No plot can be produced.")

par(mai = c(1.0, 0.75, 1.0, 0.5))

if(x$scaling == 1) {

   # Distance biplot, scaling type = 1: plot F for objects, U for variables
   # This projection preserves the Euclidean distances among the objects
   
   biplot(x$F,x$U,col=c(x$color.obj,x$color.var),xlab="PCA axis 1",ylab="PCA axis 2")
   title(main = c("PCA biplot","scaling type 1"), family="serif", line=4)

   } else {

   # Correlation biplot, scaling type = 2: plot G for objects, U2 for variables
   # This projection preserves the correlation among the variables
   
   biplot(x$G,x$U2,col=c(x$color.obj,x$color.var),xlab="PCA axis 1",ylab="PCA axis 2")
   title(main = c("PCA biplot","scaling type 2"), family="serif", line=4)

   }
invisible()
}
