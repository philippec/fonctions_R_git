CCA <- function (Y, X, use.svd=TRUE)
#
# Compute canonical correspondence analysis (CA).
# Data table Y must contain frequencies or equivalent; no negative values.
#
# Licence: GPL-2
# Author:: Pierre Legendre, November 2010
{
### Internal functions
###
stand.X <- function(X, n, m, fi., f.)
# Standardize X using the vectors of means and standard deviations of X.infl
# computed without actually constructing the inflated matrix X.infl.
# Numerical ecology (1998), p. 595, parag. 3
{
D.fi <- diag(fi.)
X.mean <- apply(D.fi %*% X, 2, sum)/f.  # Vector of means
#
temp <- sweep(X, 2, X.mean, "-")        # Matrix X centred
tmp  <- D.fi %*% temp^2                 # Matrix of squared values, weighted
X.sd <- sqrt(apply(tmp, 2, sum)/f.)     # Vector of standard deviations
#
X.stand <- sweep(temp, 2, X.sd,"/")     # Standardized matrix X
X.stand
# out <- list(X.stand=X.stand, X.mean=X.mean, X.sd=X.sd)
# out
}
###
stand.X2 <- function(X, n, m, fi., f.)
# Alternative computation method:
# Construct inflated matrix X.infl and use its vectors 
# of means and standard deviations to standardize X.
# Numerical ecology (1998), p. 595, parag. 3
{
X.infl <- matrix(NA, f., m)
no <- 0
for(i1 in 1:n) {
	for(i2 in 1:fi.[i1]) {
		no <- no+1
		X.infl[no,] <- X[i1,]
		}
	}
# Standardize X.infl
X.mean <- apply(X.infl,2,mean)
X.sd <- sqrt(apply(X.infl,2,var)*(f.-1)/f.)
#
temp <- sweep(X,2,X.mean,"-")
X.stand <- sweep(temp,2,X.sd,"/")
X.stand
# out <- list(X.stand=X.stand, X.mean=X.mean, X.sd=X.sd)
# out
}
###
w.regress <- function(Y, X, diag.w, fit=TRUE)
# Weighted regression.
# Returns Yhat (fitted values) if fit=TRUE or residuals if fit=FALSE
{
Yhat <- diag.w^(0.5) %*% X %*% solve(t(X) %*% diag.w %*% X) %*% t(X) %*% diag.w^(0.5) %*% Y
if(fit) { Yhat } else { Y-Yhat }
}
###
### End internal functions
#
Y <- as.matrix(Y)
X <- as.matrix(X)
if(min(Y) < 0) stop("CCA does not allow negative values in Y")
if(det(cor(X))==0) cat("X contains collinear variables")
if(ncol(X) > (nrow(X)-1)) stop("X contains more columns than (n-1)")
vif <- vif.JO(X)
keep <- which( !is.na(vif) )
#
# Calculate three basic parameters
n <- nrow(Y)
p <- ncol(Y)
m <- ncol(X)
#
# Save the row and column names
site.names <- rownames(Y)
sp.names <- colnames(Y)
env.names <- colnames(X)
names(vif) <- env.names
#
# Construct the Qbar matrix (contributions to chi-square)
# Numerical ecology (1998), equation 9.32
fi. <- apply(Y,1,sum)     # Row weights, f(i+)
f.j <- apply(Y,2,sum)
f.  <- sum(fi.)           # f(++)
pi. <- fi./f.
p.j <- f.j/f.
E <- ( matrix(fi.,n,1) %*% matrix(f.j,1,p) )/f.
Qbar <- (Y - E) * E^(-0.5) / sqrt(f.)
inertia <- sum(Qbar^2)
D.pi. <- diag(pi.)
# 
XX.stand <- stand.X(X, n, m, fi., f.)
X.stand <- XX.stand[,keep]
#
Yhat <- w.regress(Qbar, X.stand, D.pi., fit=TRUE)
#
if(!use.svd) {
### Analyse Yhat by eigen()
eig <- eigen(t(Yhat) %*% Yhat)
k <- length(which(eig$values > 1e-10))
values <- eig$values[1:k]
U <- eig$vectors[,1:k]
} else {
### Alternative code, svd()
svd.res <- svd(t(Yhat) %*% Yhat)
k <- length(which(svd.res$d > 1e-10))
values <- svd.res$d[1:k]
U <- svd.res$v[,1:k]
}
### Alternative code, svd()
# svd.Yhat <- svd(Yhat)
# k <- length(which(svd.Yhat$d > 1e-10))
# values <- svd.Yhat$d[1:k]^2
# U <- svd.Yhat$v[,1:k]
#
rel.values <- values / inertia
cum.rel <- cumsum(rel.values)
# Compute derived tables for triplots
Uhat <- Qbar %*% U %*% diag(values^(-0.5))
V <- diag(p.j^(-0.5)) %*% U
Vhat <- diag(pi.^(-0.5)) %*% Uhat
F <- Vhat %*% diag(values^(0.5))
Fhat <- V %*% diag(values^(0.5))
Z1 <- diag(pi.^(-0.5)) %*% Yhat %*% U
Z2 <- Z1 %*% diag(values^(-0.5))
Z.stand <- stand.X(Z1, n, k, fi., f.)  # Identical result with Z2
#
# Biplot scores, scalings 1 and 2: compute correlations of X variables with axes
biplot.sc2 <- t(XX.stand) %*% diag(pi.) %*% Z.stand
biplot.sc1 <- biplot.sc2 %*% diag(values^(0.5))
#
# Inter-set correlations of environmental variables with axes
r.spec.env <- diag(t(stand.X(F,n,k,fi.,f.)) %*% D.pi. %*% Z.stand)
cor.env.axes <- biplot.sc2 %*% diag(r.spec.env)
#
# Scaling=3
spec3 <- V %*% diag(values^(0.25))                # Species scores
Z3    <- Z1 %*% diag(values^(-0.25))              # Site scores Z
site3 <- Vhat %*% diag(values^(0.25))             # Site scores F
biplot.sc3 <- biplot.sc2 %*% diag(values^(0.25))  # Biplot scores
#
rownames(U) <- rownames(V) <- rownames(spec3) <- rownames(Fhat) <- sp.names
rownames(Uhat) <- rownames(F) <- rownames(Vhat) <- rownames(site3) <- rownames(Z1) <- rownames(Z2) <- rownames(Z3) <- site.names
ax.names <- paste("Axis",1:k,sep="")
colnames(U) <- colnames(Uhat) <- colnames(V) <- colnames(spec3) <- colnames(Vhat) <- colnames(site3) <- colnames(F) <- colnames(Fhat) <- colnames(Z1) <- colnames(Z2) <- colnames(Z3) <- colnames(cor.env.axes) <- colnames(biplot.sc1) <- colnames(biplot.sc2) <- colnames(biplot.sc3) <- ax.names
#
general <- list(inertia=inertia, values=values, rel.values=rel.values, cum.rel=cum.rel, r.spec.env=r.spec.env, cor.env.axes=cor.env.axes, vif=vif, k=k)
scaling1 <- list(species=V, sites.Z=Z1, sites.F=F, biplot=biplot.sc1)
scaling2 <- list(species=Fhat, sites.Z=Z2, sites.F=Vhat, biplot=biplot.sc2)
scaling3 <- list(species=spec3, sites.Z=Z3, sites.F=site3, biplot=biplot.sc3)
other <- list(U=U, Uhat=Uhat, X.stand=X.stand, Qbar=Qbar, site.names=site.names, sp.names=sp.names, env.names=env.names, call=match.call() )
#
out <- list(general=general, scaling1=scaling1, scaling2=scaling2, scaling3=scaling3, other=other)
class(out) <- "CCA"
out
}

`vif.JO` <-       # J. Oksanen function vif.cca, modified
    function(XX) 
{
    XX <- scale(XX, center=TRUE, scale=FALSE)
    Q <- qr(XX)
    out <- rep(NA, NCOL(Q$qr))
    names(out)[Q$pivot] <- colnames(Q$qr)
    rank <- Q$rank
    V <- chol2inv(Q$qr, size = rank)
    X <- qr.X(Q)[, Q$pivot[1:rank], drop=FALSE]
    Vi <- crossprod(X)
    v1 <- diag(V)
    v2 <- diag(Vi)
    out[Q$pivot[1:rank]] <- v1 * v2
    out
}

`print.CCA` <-    # Prints basic statistics of eigen-decompostion 
    function(x, ...)
{
if (!inherits(x, "CCA")) stop("Object of class 'CCA' expected")
    cat("\nCanonical Correspondence Analysis\n")
    cat("\nCall:\n")
    cat(deparse(x$other$call),'\n')
    cat("\nTotal inertia in matrix Qbar: ",x$general$inertia,'\n')
    cat("\nCanonical eigenvalues",'\n')
    cat(x$general$values,'\n')
    cat("\nRelative eigenvalues",'\n')
    cat(x$general$rel.values,'\n')
    cat("\nCumulative relative eigenvalues",'\n')
    cat(x$general$cum.rel,'\n')
    cat("\nVariance inflation factors (VIF)",'\n')
    print(x$general$vif)
    cat('\n')
    invisible(x) 
}

`triplot.CCA` <-    # Produces a triplot
    function(x, xax=1, yax=2, scaling=1, graph.type="Z", aspect=1, cex=1, color.sites="black", color.sp="red", color.env="blue", ...)
# xax and yax determine the axes that will be plotted
# Use aspect=NA to remove the effect of parameter 'asp' in the graphs
{
if (!inherits(x, "CCA")) stop("Object of class 'CCA' expected")
if(length(x$general$values) < 2) stop("There is a single eigenvalue. No plot can be produced.")
#
sp.names  <- x$other$sp.names
si.names  <- x$other$site.names
env.names <- x$other$env.names
#

if(scaling == 1) {

# The sites are at the centroids (barycentres) of the species
# This projection preserves the chi-square distance among the sites
label <- "scaling type 1"
sp <- x$scaling1$species
bi <- x$scaling1$biplot
if(graph.type == "Z") {si <- x$scaling1$sites.Z} else {si <- x$scaling1$sites.F}

} else if(scaling == 2) {

# The species are at the centroids (barycentres) of the sites
# This projection preserves the chi-square distance among the species
label <- "scaling type 2"
sp <- x$scaling2$species
bi <- x$scaling2$biplot
if(graph.type == "Z") {si <- x$scaling2$sites.Z} else {si <- x$scaling2$sites.F}

} else if(scaling == 3) {

# Biplot, scaling = 3 (Symmetric scaling in Canoco)
label <- "scaling type 3"
sp <- x$scaling3$species
bi <- x$scaling3$biplot
if(graph.type == "Z") {si <- x$scaling3$sites.Z} else {si <- x$scaling3$sites.F}

} else {

stop("Program stopped: error in scaling type")
}

# Find the limits of the axes
si.range <- apply(si[,c(xax,yax)],2,range)
sp.range <- apply(sp[,c(xax,yax)],2,range)
bi.range <- apply(bi[,c(xax,yax)],2,range)

# Plot 'sp' for species, 'si' for sites, and 'bi' for explanatory variables
ran.sp <- sp.range[2,] - sp.range[1,]
ran.si <- si.range[2,] - si.range[1,]
ran.bi <- bi.range[2,] - bi.range[1,]
ran.x <- max(ran.si[1], ran.sp[1], ran.bi[1])
fact <- min(round(ran.x/ran.bi))
xmin <- min(sp.range[1,1], si.range[1,1], bi.range[1,1]*fact) - ran.x/8
xmax <- max(sp.range[2,1], si.range[2,1], bi.range[2,1]*fact) + ran.x/3
ymin <- min(sp.range[1,2], si.range[1,2], bi.range[1,2]*fact)
ymax <- max(sp.range[2,2], si.range[2,2], bi.range[2,2]*fact)
#
plot(si[,c(xax,yax)], asp=aspect, pch=20, cex=cex, xlim=c(xmin,xmax), ylim=c(ymin,ymax), xlab=paste("CCA axis",xax), ylab=paste("CCA axis",yax), col=color.sites)
text(si[,c(xax,yax)], labels=si.names, cex=cex, pos=4, offset=0.5, col=color.sites)
points(sp[,c(xax,yax)], pch=22, cex=cex, col=color.sp)
text(sp[,c(xax,yax)], labels=sp.names, cex=cex, pos=4, offset=0.5, col=color.sp)
arrows(x0=0,y0=0, x1=bi[,xax]*fact, y1=bi[,yax]*fact, lty=1, col=color.env, code=2, len=0.1, lwd=1)
text(bi[,c(xax,yax)]*fact, labels=env.names, cex=cex, pos=4, offset=0.5, col=color.env)
title(main <- c("CCA triplot",label), family="serif")
#
invisible()
}