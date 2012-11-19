WRperiodogram <- function(x, T1=2, T2=NULL, nperm=499, mult="bonferroni")
#
# Whittaker and Robinson periodogram; Whittaker and Robinson (1924),
# Legendre & Legendre (1998, 2012, Section 12.4.1).
# x : a vector of quantitative data
# T1: first period to analyse
# T2: last period to analyse (T2 <= n/2)
# nperm: number of permutations for tests of significance
# mult : methods for correction for multiple testing; "sidak" or "bonferroni"
#
# License: GPL-2 
# Author:: Pierre Legendre, September 2012
{
### Internal function
###
BB <- function(x,T1,T2,n,nk)
# Construct Buys-Ballot tables for periods T1 to T2.
# For each table (i.e. for each period), compute the periodogram statistic 'sd', 
# which is the square-root of the variance of the column means.
# Use the maximum-likelihood estimator of the variance, with division by n,
# instead of (n-1) in the ordinary unbiased estimator of the variance.
# WR.stat is the Whittaker-Robinson standard deviation periodogram statistic.
{
	WR.stat <- NA
	for(k in T1:T2) {
		nr <- ceiling(n/k)
		vec1 <- c(x,rep(NA,(nr*k-n)))
		tab <- matrix(vec1, nr, k, byrow=TRUE)
		mean.vec <- apply(tab, 2, mean, na.rm=TRUE)
		sel <- which(is.na(mean.vec))
		if(length(sel)>0) kk <- (k-length(sel)) else kk <- k
		WR.stat <- c(WR.stat, sqrt(var(mean.vec, na.rm=TRUE)*(kk-1)/kk))
		}
	WR.stat[-1]
}
### End internal function

if(is.factor(x)) stop("x contains a factor; it sould contain real numbers")
mult <- match.arg(mult, c("sidak", "bonferroni"))
# Make all values positive: this is not necessary for this periodogram
# if(min(x)<0) x <- x - min(x)
#
n <- length(x)
T.sup <- n %/% 2
if(is.null(T2)) T2 <- T.sup
if(T2 > T.sup) T2 <- T.sup
nk <- T2-T1+1   # Number of tests performed
# cat("T1 =",T1,"T2 =",T2,"nk =",nk,"\n")
out <- matrix(NA,nk,5)
colnames(out) <-c("Period","WR.stat","p.value","p.corrected","p.corr.progr")

WR.stat <- BB(x,T1,T2,n,nk)
out[,1] <- T1:T2
out[,2] <- WR.stat

# Permutation test on the 'nk' statistics
if(!is.null(nperm) & nperm>0) {
	prob <- rep(1,nk)   # Holm (1967) correction
	for(i in 1:nperm) {
		WR.perm <- BB(sample(x),T1,T2,n,nk)
		prob <- prob + (WR.perm >= WR.stat)
		# cat("prob =",prob,"\n")
		}
		prob <- prob/(nperm+1)
	out[,3] <- prob
	# Correct all p-values for 'n.tests' simultaneous tests
	    if(mult == "sidak") {
	    	# Sidak correction
	        p.corr <- 1 - (1 - prob)^nk
	        # Sidak with progressive p-value correction
	        p.corr.prog <- 1 - (1 - prob)^(1:nk)
	    } else { 
	    	# Bonferroni correction
	        p.corr <- prob*nk
	        # Bonferroni with progressive p-value correction
	        p.corr.prog <- prob*(1:nk)
	    }
	p.corr[which(p.corr>1)] <- 1
	out[,4] <- p.corr
	p.corr.prog[which(p.corr.prog>1)] <- 1
	out[,5] <- p.corr.prog
	}
# Check for NA in output file
if(any(is.na(out[,2]))) {
	sel <- which(is.na(out[,2]))
	sel.true <- TRUE
	cat("The following period(s) removed from output file due to NA:",sel+1,"\n")
	} else { 
	sel.true <- FALSE
	}
if(nperm>0) {
	sel.na <- which(is.na(out[,3]))
	out[sel.na,3:5] <- 99
	}
if(sel.true) out <- out[-sel,]
class(out) <- "WRperio"
out
}

plot.WRperio <- function(x, prog=1, alpha=0.05, line.col="red")
# Plot a Whittaker-Robinson periodogram
# prog = 1: use the orginal p-values in the plot.
#      = 2: use the p-values corrected for multiple testing.
#      = 3: progressive correction of multiple-testing (default).
# alpha: significance level for the plot. Default: alpha=0.05.
{
plot(x[,1], x[,2], col="red", xlab="Period", 
	ylab="Whittaker-Robinson stat.", ylim=c(0, max(x[,2])))
lines(x[,1], x[,2], col=line.col)              # Plot W-R statistics
points(x[,1], x[,2], pch=22, bg="white")       # Re-plot points as white squares
if(is.na(x[1,3])) cat("p-values in the W-R periodogram were not computed\n")
select <- which(x[,(prog+2)] <= alpha)
points(x[select,1], x[select,2], pch=22, bg="black") # Significant points: black
}
