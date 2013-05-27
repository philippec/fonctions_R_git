'cross' <- 
function(Y, res, k1=2, k2=NULL, xv="min", xvmult=100)
# Construct a file with the membership of all clustering levels (columns).
# That file will serve as the constraining data in mvpart().
#
# res = the constrained.clust or hclust output object.
# k1 = Lowest number of clusters (default: k1=2) for group coding.
# k2 = Highest number of clusters for group coding. If k2=NULL, k2=nlev (highly recommended since the best solution may involve small groups, which are found in classifications with high k2 values).
#
# License: GPL-2 
# Author:: Pierre Legendre, January 2011
{
### Internal function
"simpleRDA" <- function (Y, X, SS.Y, ...)
{
    Q <- qr(X, tol=1e-6)
    Yfit.X <- qr.fitted(Q, Y)
    SS <- sum(Yfit.X^2)
    if (missing(SS.Y)) SS.Y <- sum(Y^2)
    Rsquare <- SS/SS.Y
    list(Rsquare = Rsquare, m = Q$rank)
}
### End internal function
nlev <- nrow(res$merge)
if( is.numeric(Y) & !(is.matrix(Y)) ) Y = cbind(Y, rep(1,length(Y)))
if( (is.matrix(Y) | is.data.frame(Y)) & ncol(Y)==1 ) 
	Y = cbind(Y[,1], rep(1,nrow(Y)))
Y <- scale(Y, center=TRUE, scale=FALSE)
n <- nrow(Y)
if( is.null(k2) ) {
	k2 <- ceiling(nlev/2)
	} else { 
	if(k2>nlev) {
		k2 <- ceiling(nlev/2)
		cat("k2 > No. levels in hierarchical classification","\n")
		cat("Corrected value of k2 =",k2,"\n")
		}
	}
if(k2 > 15) {
	cat("Are you sure you want cross-validation for",k2,"groups?\n")
	cat("Calculations may take a long time.\n")
	cat("Type 'Y' if yes, 'N' if no.\n")
	answer <- toupper(scan(file="",what="character",nlines=1,quiet=TRUE))
	if((answer=="n") | (answer=="N")) stop("Calculation interrupted by user")
	}
#
# Construct data.frame "out" containing factors representing the partitions
out <- as.data.frame(as.factor(cutree(res, k1)))
for(j in (k1+1):k2) { out <- cbind(out, as.factor(cutree(res, j))) }
# cat("k1 =", k1, " k2 =", k2, "\n")
# cat("nrow(out) =", nrow(out), " ncol(out) =", ncol(out), "\n")
if(length(rownames(Y)) == 0) {
	rownames(out) = rownames(Y)
	} else {
	rownames(out) <- paste("obj.", 1:n, sep="")
	}
colnames(out) <- paste("Gr.", k1:k2, sep="")

# Compute AICc and Calinski-Harabasz criterion for each partition
AIC <- matrix(NA, (k2-k1+1), 6)
rownames(AIC) <- colnames(out)
colnames(AIC) <- c("Rsquare", "AICc", "C-H", "P(C-H)", "cvre", "ngr.cross")

for(j in 1:(k2-k1+1)) {
	# cat("Partition column #",(k1-no.1+j),"\n")
	target <- as.data.frame(as.factor(out[,j]))
	target <- model.matrix(~.,data=target)  # Keeping the intercept = centring
	temp <- simpleRDA(Y, target)
	AIC[j,1] <- R2 <- temp$Rsquare
	m <- ncol(target)-1

#	temp <- simpleRDA(Y, out[,j])
#	R2 <- temp$Rsquare
#	m <- length(levels(out[,j])) - 1
	# cat("Partition column #",(k1-no.1+j)," n =",n," m =",m,"\n")
	AIC[j,2] <- log((1-R2)/n) + (n+m)/(n-m-2)  # AICc
	AIC[j,3] <- R2*(n-1-m) / ((1-R2)*m)        # Calinski-Harabasz criterion 
	AIC[j,4] <- pf(AIC[j,3], m, (n-1-m), lower.tail=FALSE)  # P of C-H F-stat.
	}

# Cross-validation using {mvpart} on each partition in 'out', one at a time
for(j in 1:(k2-k1+1)) {
	cat("*** Cross validation, solution with",(k1+j-1),"groups",'\n')
	invisible(capture.output( MRT.res <- mvpart(data.matrix(Y) ~ out[,j], 
		data=as.data.frame(out[,j]), xv=xv, xvmult=xvmult) ))
	cptable.res <- MRT.res$cptable
	AIC[j,6] <- ngr <- nrow(MRT.res$cptable)
	AIC[j,5] <- cptable.res[ngr,4]
	}
	cat("*** End of cross validation",'\n')

# Find lowest value of AICc, P and cvre
AIC.min <- which( AIC[,2] == min(AIC[which(!is.nan(AIC[,1])),2]) )
P.min <- which( AIC[,4] == min(AIC[which(!is.nan(AIC[,4])),4]) )
cvre.min <- which( AIC[,5] == min(AIC[which(!is.nan(AIC[,5])),5]) )
cvre.col <- AIC[cvre.min,6]-k1+1

#
res <- list(AIC=AIC, ngr.AIC=(AIC.min+k1-1), classif.AIC=out[,AIC.min], ngr.P=(P.min+k1-1), classif.P=out[,P.min], cvre.and.ngr=c(AIC[cvre.min,5:6]), classif.cvre=out[,cvre.col], out=out)
res
}
