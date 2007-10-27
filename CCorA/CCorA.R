CCorA <- function(Y, X1, X2=NULL, stand.Y=FALSE, stand.X1=FALSE, stand.X2=FALSE,
         print.plot=TRUE, print.obj=FALSE)
#
################################################################################
#
# Canonical and partial canonical correlation analysis, Y <-> X1/X2. 
# See the help file "CCorA.help.pdf".
#
# This function will be incorporated into the vegan library.
#
#                                                     Pierre Legendre, June 2006
#
################################################################################
{
library(MASS)
if(is.logical(stand.Y)){
	} else { stop("Wrong operator; 'center.Y' should be either 'FALSE' or 'TRUE'") }
if(is.logical(stand.X1)){
	} else { stop("Wrong operator; 'center.X1' should be either 'FALSE' or 'TRUE'") }
if(is.logical(stand.X2)){
	} else { stop("Wrong operator; 'center.X2' should be either 'FALSE' or 'TRUE'") }
partial <- FALSE

Y <- as.matrix(Y)
var.null(Y,1)
nY <- nrow(Y)
p <- ncol(Y)
X1 <- as.matrix(X1)
var.null(X1,2)
nX1 <- nrow(X1)
q <- ncol(X1)
if(nY != nX1) stop("Program stopped: Different numbers of rows in Y and X1")
n <- nY
if((p+q) >= (n-1)) stop("Program stopped: Not enough degrees of freedom!")
rownoms <- rownames(Y)

Y.c <- apply(Y,2,scale,center=TRUE,scale=stand.Y)
X1.c <- apply(X1,2,scale,center=TRUE,scale=stand.X1)

# Replace Y.c and X1.c by tables of their PCA object scores, computed by SVD
temp <- cov.inv(Y.c,1)
Y <- temp$mat
S.Y.inv <- temp$S.inv
# colnames(Y) <- colnames(Y.c[,1:temp$m])
rownames(Y) <- rownoms
if((temp$m == 1) & (print.plot == TRUE)) {
   print.plot <- FALSE
   cat("No plot will be produced because Y has a single dimension",'\n')
   }
#
temp <- cov.inv(X1.c,2)
X1 <- temp$mat
S.X1.inv <- temp$S.inv

if(length(X2) != 0) {   # Compute residuals of X1 over X2
   partial <- TRUE
   X2 <- as.matrix(X2)
   var.null(X2,3)
   X2.c <- apply(X2,2,scale,center=TRUE,scale=stand.X2)
   nX2 <- nrow(X2.c)
   if(nY != nX2) stop("Program stopped: Different numbers of rows in Y and X2")
   # Replace X2.c by the table of its PCA object scores, computed by SVD
   temp2 <- cov.inv(X2.c,3)
   X2 <- temp2$mat
   }

if(partial == FALSE) {
   X <- X1
   S.X.inv <- S.X1.inv
   if((temp$m == 1) & (print.plot == TRUE)) {
      print.plot <- FALSE
      cat("No plot will be produced because X1 has a single dimension",'\n')
      }
   } else {
   # Regress X1 on X2 and compute res(X1)/X2 before the canonical analysis to estimate [a]
   Q <- qr(X2, tol=1e-6)
   X <- qr.resid(Q, X1)
   # Replace X by the table of its PCA object scores, computed by SVD
   temp <- cov.inv(X,4)
   X <- temp$mat
   S.X.inv <- temp$S.inv
   if((temp$m == 1) & (print.plot == TRUE)) {
      print.plot <- FALSE
      cat("No plot will be produced because X1 has a single dimension after controlling for X2",'\n')
      }
   }
# colnames(X) <- colnames(X1.c[,1:temp$m])
rownames(X) <- rownoms

# Covariance matrices, etc. from the PCA scores
epsilon <- 0.0000000001
S11 <- cov(Y)
if(sum(abs(S11)) < epsilon) return(0)
S22 <- cov(X)
if(sum(abs(S22)) < epsilon) return(0)
S12 <- cov(Y,X)
if(sum(abs(S12)) < epsilon) return(0)

S11.chol <- chol(S11)
S11.chol.inv <- solve(S11.chol)
S22.chol <- chol(S22)
S22.chol.inv <- solve(S22.chol)

# K summarizes the correlation structure between the two sets of variables
K <- t(S11.chol.inv) %*% S12 %*% S22.chol.inv

K.svd <- svd(K)
EigenValues <- K.svd$d^2
# K.svd$u %*% diag(K.svd$d) %*% t(K.svd$v)   # This line checks that K = U D V'
axenames <- paste("CanAxis",1:length(K.svd$d),sep="")
U <- K.svd$u
V <- K.svd$v

A <- S11.chol.inv %*% U
B <- S22.chol.inv %*% V

Cy <- (Y %*% A)/sqrt(n-1)
Cx <- (X %*% B)/sqrt(n-1)

# Compute the 'Biplot scores of Y variables' a posteriori --
# use 'ginv' for inversion in case there is collinearity
# AA <- coefficients of the regression of Cy on Y.c, times sqrt(n-1)
# AA <- sqrt(n-1) * [Y'Y]-1 Y' Cy
YprY <- t(Y.c) %*% Y.c
AA <- sqrt(n-1) * ginv(YprY) %*% t(Y.c) %*% Cy
#
# Compute the 'Biplot scores of X variables' a posteriori --
XprX <- t(X1.c) %*% X1.c
BB <- sqrt(n-1) * ginv(XprX) %*% t(X1.c) %*% Cx
	
rownames(U) <- colnames(Y)
rownames(V) <- colnames(X)
rownames(Cy) <- rownoms
rownames(Cx) <- rownoms
colnames(U) <- axenames
colnames(A) <- axenames
colnames(AA) <- axenames
colnames(V) <- axenames
colnames(B) <- axenames
colnames(BB) <- axenames
colnames(Cy) <- axenames
colnames(Cx) <- axenames
# pp <- length(K.svd$d)

# Check U and V by eigenvalue decomposition
# KKpr.eig <- eigen(K %*% t(K))
# eigval1 <- KKpr.eig$values
# UU <- KKpr.eig$vectors
# KprK.eig <- eigen(t(K) %*% K)
# eigval2 <- KprK.eig$values
# VV <- KprK.eig$vectors

# Compute Pillai's trace = sum of the canonical eigenvalues
#                        = sum of the squared canonical correlations
gross.mat <- S12 %*% solve(S22) %*% t(S12) %*% solve(S11)
PillaiTrace <- sum(diag(gross.mat))

# Graphs
if(print.plot == TRUE) {
   par(mfrow=c(1,2))
   if(print.obj == TRUE) {
      biplot(Cy, AA)
      biplot(Cx, BB)
      } else {
      biplot(Cy, AA, col=c("white","black"))
      biplot(Cx, BB, col=c("white","black"))
      }      
   }
# Compute the two redundancy statistics
RsquareY.X <- simpleRDA2(Y, X)
RsquareX.Y <- simpleRDA2(X, Y)
Rsquare.adj.Y.X <- RsquareAdj(RsquareY.X$Rsquare, n, RsquareY.X$m)
Rsquare.adj.X.Y <- RsquareAdj(RsquareX.Y$Rsquare, n, RsquareX.Y$m)

return(list(Pillai=PillaiTrace, EigenValues=EigenValues, CanCorr=K.svd$d,
   Mat.ranks=c(RsquareX.Y$m, RsquareY.X$m), 
   RDA.Rsquares=c(RsquareY.X$Rsquare, RsquareX.Y$Rsquare),
   RDA.adj.Rsq=c(Rsquare.adj.Y.X, Rsquare.adj.X.Y),
#      eigval1=eigval1, eigval2=eigval2, UU=UU, VV=VV,
#      U=U, A=A, V=V, B=B,
   AA=AA, BB=BB, Cy=Cy, Cx=Cx))
}

simpleRDA2 <- function (Y, X, SS.Y, ...)  # This function is already in vegan
{
    Q <- qr(X, tol = 1e-06)
    Yfit.X <- qr.fitted(Q, Y)
    SS <- sum(Yfit.X^2)
    if (missing(SS.Y)) 
        SS.Y <- sum(Y^2)
    Rsquare <- SS/SS.Y
    list(Rsquare = Rsquare, m = Q$rank)
}

RsquareAdj <- function (x, n, m, ...) # This function is already in vegan
{
    if (m >= (n - 1)) 
        NA
    else 1 - (1 - x) * (n - 1)/(n - m - 1)
}

var.null <- function(mat, no)
# Check for the presence of variables with null variances
{
mat.cov <- cov(mat)
problem <- FALSE
for(i in 1:nrow(mat.cov)) {
   if(mat.cov[i,i] == 0) {
      cat("Matrix",no,"-- Variable no.",i," has a null variance",'\n')
      problem <- TRUE
      }
   }
if(problem == TRUE) stop("Program stopped. Verify/modify your matrix No.",no)
}

cov.inv <- function(mat, no)
#
# This function returns:
#
# 1) mat = the matrix of PCA object scores (by SVD);
# 2) S.inv = the inverse of the covariance matrix;
# 3) m = the rank of matrix 'mat'
#
# The inverse of the PCA covariance matrix is simply the diagonal matrix of (1/eigenvalues).
# If ncol(mat) = 1, the inverse of the covariance matrix simply contains 1/var(mat).
{
mat <- as.matrix(mat)
if(ncol(mat) == 1) {
   S.inv <- as.matrix(1/var(mat))
   m <- 1
   } else {
   epsilon <- 0.00000001
   S.svd <- svd(cov(mat))
   m <- ncol(mat)
   mm <- 0
   for(i in 1:m) { if(S.svd$d[i] > epsilon) mm <- mm+1 }
   if(mm < m) {
      cat("Matrix",no,"  S: rank=",mm," < order",m,'\n')
      if((mm == 0) & (no == 4)) stop("Program stopped: X1 has rank = 0 after controlling for X2")
      m <- mm
      }
   S.inv <- diag(1/S.svd$d[1:m])
   mat <- mat %*% S.svd$u[,1:m]
   }
return(list(mat=mat, S.inv=S.inv, m=m))
}
