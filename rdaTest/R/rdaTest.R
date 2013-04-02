'rdaTest' <- 
	function(Y, X, W=NULL, scale.Y=FALSE, test.F=NULL, nperm=NULL, algorithm=1, 
	silent=FALSE, seed=NULL, verbose=FALSE)
#
# This function computes partial canonical redundancy analysis (partial RDA).
#
# Parameters --
#
#    algorithm=1: partial RDA is RDA of Y with respect to X(res|W) 
#    algorithm=2: partial RDA is RDA of Y(res|W) with respect to X(res|W) 
#
#    seed: seed for random number generator, used by the permutation function 
#          sample(). If seed=NULL (default), a random integer is drawn as the 
#          seed for the run. It will be reset to that value before the test of 
#          each canonical axis. All axes are thus tested using the same set of 
#          permutations. A fixed value, e.g. seed=12345, can be givenn by the 
#          user to compare the results of this function with that of other 
#          functions where the seed can also be set at run time.
#    verbose = TRUE: print intermediate F and F.perm results
#
# License: GPL-2 
# Authors: Pierre Legendre, using ideas of C.J.F. ter Braak and J. Oksanen, 2010
{
	if(length(seed)==0) seed <- ceiling(runif(1,max=10000))
	cat("seed =",seed,'\n')

	if(!is.logical(scale.Y)) 
		stop("Wrong operator; 'scale.Y' should be either 'FALSE' or 'TRUE'")	
	epsilon <- sqrt(.Machine$double.eps)	
	# Read the data tables. Transform them into matrices Y, X, and W
	Y.mat = as.matrix(Y)
	X.mat = as.matrix(X)
	sitenames <- rownames(X.mat)
	if(is.null(rownames(X.mat)) & is.null(rownames(Y.mat))) {
		sitenames <- paste("Site",1:nrow(X.mat),sep="")
		} else {
			if(is.null(rownames(X.mat))) {
				sitenames <- rownames(Y.mat)
			} else {
				sitenames <- rownames(X.mat)
			}
		}
	
	# If(scale.Y == TRUE), check the presence of Y.mat variables 
	# with null variances
	if(scale.Y) {
		mat.sd <- apply(Y.mat,2,sd)
		problem = FALSE
		for(i in 1:length(mat.sd)) {
		   if(mat.sd[i] == 0) {
		      cat(" Variable No.",i," in Y has a null variance",'\n')
		      problem = TRUE
		      }
		   }
		if(problem == TRUE) 
			stop("The program was stopped. Verify/modify the Y matrix.")	
		}

	# Check the presence of X.mat variables with null variances
	mat.sd <- apply(X.mat,2,sd)
	problem = FALSE
	for(i in 1:length(mat.sd)) {
	   if(mat.sd[i] == 0) {
	      cat(" Variable No.",i," in X has a null variance",'\n')
	      problem = TRUE
	      }
	   }
	if(problem == TRUE) 
		stop("The program was stopped. Verify/modify the X matrix.")
	
	# Pre-treatment: centre or standardize Y, standardize X
	n = nrow(Y.mat)
	p = ncol(Y.mat)
	Y = scale(Y.mat,center=TRUE,scale=scale.Y)    ### Y
	YY = Y                                        ### YY
	rownames(Y) = sitenames
	SS.Y <- sum(Y^2)                          # [a+b+c+d]
	X = scale(X.mat,center=TRUE,scale=TRUE)
	XX = X
#	cat(colnames(X),'\n')
	
	# Identify collinear columns in X and remove them
	m1 = ncol(X)
	if(m1 > 1) {
		vif.res = vif(X, print.vif=!silent)
		if(length(which(vif.res==0))!=0) X = X[,-which(vif.res==0)]
		} else { 
		vif.res = rep(1,m1) 
		}
	m = ncol(as.matrix(X))
	
	if((!silent) & m > 1) { 
		cat('\n')
		cat('Variance Inflation Factors (VIF)','\n')
		print(vif.res)
	}		
#	cat(colnames(X),'\n')

	qq = 0
	mq = m
	# Is there a matrix W containing covariables?
	if(!is.null(W)) {
		covar = TRUE
		W.mat = as.matrix(W)
		qq = ncol(W.mat)
		W = scale(W.mat,center=TRUE,scale=TRUE)
    	if(!silent) cat("\nThere is a matrix of covariables\n")
    	# Find the rank of W 
    	q <- sum(svd(cov(W), nv = 0, nu = 0)$d > epsilon)

        # Find the rank of cbind(X,W) 
        mq <- sum(svd(cov(cbind(X, W)), nv = 0, nu = 0)$d > epsilon)
		} else {
		covar=FALSE
    	if(!silent) cat("\nThere is no matrix of covariables\n")
		}

	# If covariables W are present, compute Y.res.W and X.res.W
	if(covar) {
        qr.W <- qr(W)
        Y <- qr.resid(qr.W, YY)        # Y(res|W)       ### Y = [a+d]
		SS.Yres.W = sum(Y^2)           # [a+d]
		X <- qr.resid(qr.W, X)         # X(res|W)
	
		XW = cbind(X,W)
		qr.XW <- qr(XW)
		XW.fit <- qr.fitted(qr.XW, YY)  # We need to use YY here

		SS.Yfit.XW = sum(XW.fit*XW.fit) # [a+b+c]
		}
	
	qr.X <- qr(X)
	Yfit.X <- qr.fitted(qr.X, Y)
	rownames(Yfit.X) <- sitenames
	
	# Compute R-square and adjusted R-square
	SS.Yfit.X = sum(Yfit.X^2)           # [a]
	Rsquare = SS.Yfit.X / SS.Y
	if(!covar) {
		totalDF = n-1
		residualDF = n-m-1
		if(residualDF > 0) {
			adjRsq = 1-((1-Rsquare)*totalDF/residualDF)
			} else {
			adjRsq = NA
			if(!silent) cat("\nAdjusted R-square not computed: ",
			"residual d.f. <= 0\n")
			}
		}

	# Test significance of the canonical F statistic
	testFF="N"
	if(is.null(test.F)) {
		cat('\n')
		cat("Test the F-statistic? Type 'Y' if yes, 'N' if no.",'\n')
		testFF=toupper(scan(file="",what="character",nlines=1,quiet=TRUE))
	    if(testFF=="y") testFF="Y"
		} else { 
		if(test.F==TRUE) testFF="Y" 
		}

	if(testFF=="Y") {
		if(is.null(nperm)) {
			cat('\n','How many permutations? Ex. 499, 999, ...','\n',sep="")
			nper=toupper(scan(file="",what="integer",nlines=1,quiet=TRUE))
			nperm=as.integer(nper)
			}
		if(nperm > 5) verbose <- FALSE   # Modify this limit as needed #
		}
	if(!silent) {
   		cat('\n','---','\n','\n',sep="")
		cat('Redundancy statistic (canonical R-square):','\n')
		if(covar) {
			cat('Semipartial R-square =',Rsquare,
			'-- No adjusted R-square in partial RDA','\n')
			} else {
			cat('R-square =',Rsquare)
			cat(';   adjusted R-square = ',adjRsq,'\n')
			}
		}

	a <- system.time({             # How much time for the permutation test?
	if(testFF=="Y") {
	prob <- probFrda(YY,X,n,p,m,mq,nperm,qr.X,qr.W,qr.XW,SS.Y,SS.Yfit.X,
	    SS.Yfit.XW,covar, silent=silent, seed=seed, verbose=verbose)
		} else { 
		prob = "Not tested" 
		}
	})
	a[3] <- sprintf("%2f",a[3])
	if(!silent) cat('\nTime for permutation test =',a[3]," sec",'\n')
	
	# PCA portion of RDA: eigenanalysis, then F.sc1 and Z.sc1 matrices, etc.
	Y.fit.eig <- svd(Yfit.X)            # svd decomposition of Y.fit
	eig.values <- (Y.fit.eig$d^2)       # Same as eigen(cov(Y.fit))$values*(n-1)
	k <- sum(eig.values > epsilon)      # Number of canonical eigenvalues
	eig.values <- eig.values[1:k]       
	axenames <- paste("Axis",1:k,sep="")
	names(eig.values) <- axenames
	EigvalPercent = 100*eig.values/SS.Y
	cumulVar = cumsum(EigvalPercent)

	# For (scaling==1)
	U <- as.matrix(Y.fit.eig$v[,1:k])     # Same as eigen(cov(Y.fit))$vectors
	F.sc1 = Y %*% U
	Z.sc1 = Yfit.X %*% U
	
	# For (scaling==2)
	if(k == 1) {
		U.sc2 <- U * sqrt(eig.values[1]/(n-1))
   		F.sc2 = F.sc1 * (1/sqrt(eig.values[1]/(n-1)))
   		Z.sc2 = Z.sc1 * (1/sqrt(eig.values[1]/(n-1)))
   		} else {
		U.sc2 <- U %*% diag(sqrt(eig.values/(n-1)))
		F.sc2 = F.sc1 %*% diag(1/sqrt(eig.values/(n-1)))
		Z.sc2 = Z.sc1 %*% diag(1/sqrt(eig.values/(n-1)))
   		}
	var.names <- colnames(Y,do.NULL=FALSE,prefix="Var")
	rownames(U)     <- rownames(U.sc2) <- var.names
	rownames(F.sc1) <- rownames(Z.sc1) <- sitenames
	rownames(F.sc2) <- rownames(Z.sc2) <- sitenames
	colnames(U)     <- colnames(F.sc1) <- colnames(Z.sc1) <- axenames
	colnames(U.sc2) <- colnames(F.sc2) <- colnames(Z.sc2) <- axenames

	# Compute the 'Biplot scores of environmental variables' --
	#
	# Biplot scores for scaling 1: 
	# First, compute the correlations between X and Z.sc1, or X and Z.sc2 ...

	corXZ = cor(XX, Z.sc1)
	
	# ...then, weigh these correlations by the diagonal matrix 'D' of weights
	# Scale the arrows using SS.Y or, if(covar), SS.Yres.W
	# The arrow directions are the same with both scaling constants
	if(k == 1) { 
   	   if(covar) {
   	      posX.sc1 = corXZ * (sqrt(eig.values[1]/SS.Yres.W))
   	      } else {
   	      posX.sc1 = corXZ * (sqrt(eig.values[1]/SS.Y))
   	      }
	   } else {
   	   if(covar) {
   	      D = diag(sqrt(eig.values/SS.Yres.W))
   	      } else {
   	      D = diag(sqrt(eig.values/SS.Y))
   	      }
   	   posX.sc1 = corXZ %*% D
   	   }
   	
	# Biplot scores for scaling 2:
	posX.sc2 = corXZ                      # = cor(X,Z.sc1) = cor(X,Z.sc2)
	
	# If k = 1, compute the first PCA of the residuals and add it to all tables
	if(k == 1) {
	   Yres.X = Y - Yfit.X                # As in rda() {vegan}
       Yres.eig <- svd(Yres.X)            # svd decomposition of Yres.X
       Yres.values <- (Yres.eig$d^2)      # = eigen(cov(Yres.X))$values*(n-1)
       k1 <- sum(Yres.values > epsilon)

       Yres.U <- as.matrix(Yres.eig$v[,1:k1])        # PCA eigenvectors 
       Yres.F <- as.matrix(Yres.X %*% Yres.U)

       Yres.U2 <- as.matrix(Yres.U[,1]) * (Yres.values[1]^(0.5))
       Yres.G  <- as.matrix(Yres.F[,1]) * (Yres.values[1]^(-0.5))

       eig.values = c(eig.values, Yres.values[1])
       U = cbind(U, Yres.eig$v[,1])
       F.sc1 = cbind(F.sc1, Yres.F[,1])
       Z.sc1 = cbind(Z.sc1, Yres.F[,1])
       U.sc2 = cbind(U.sc2, Yres.U2)
       F.sc2 = cbind(F.sc2, Yres.G)
       Z.sc2 = cbind(Z.sc2, Yres.G)
       posX.sc1 = cbind(posX.sc1, 0)
       posX.sc2 = cbind(posX.sc2, 0)
       colnames(U) <- colnames(U.sc2) <- c("Axis1","PCA_Axis1")
       colnames(F.sc1) <- colnames(F.sc2) <- c("Axis1","PCA_Axis1")
       colnames(Z.sc1) <- colnames(Z.sc2) <- c("Axis1","PCA_Axis1")
       colnames(posX.sc1) <- colnames(posX.sc2) <- c("Axis1","PCA_Axis1")
	   }
	if(m==1) {
		rownames(posX.sc1) <- rownames(posX.sc2) <- "Expl"
		} else {
		rownames(posX.sc1) <- colnames(XX,do.NULL=FALSE,prefix="Expl")
		rownames(posX.sc2) <- colnames(XX,do.NULL=FALSE,prefix="Expl")
		}
	if(k>1) colnames(posX.sc1) <- colnames(posX.sc2) <- axenames
	
	# Print results
   	if(!silent)	{
   		cat('\n','Number of objects:  n = ',n,'\n',sep="")
   		cat('Number of response variables in Y:  p =',p,'\n')
   		cat('Number of explanatory variables in X: ',m1,
   			';  Rank of X:  m =',m,'\n')
   		if(covar==TRUE) {
   		   cat('Number of covariables in W: ',qq,';  Rank of W:  q =',q,'\n')
   		   cat('Rank of cbind(X,W):  mq =',mq,'\n')
   		   }
   		cat('Number of canonical eigenvalues:  k =',k,'\n')
   		cat('\n')
		cat('Total variance =',SS.Y/(n-1),'\n','\n')
		cat('Eigenvalues','\n')
		cat(eig.values/(n-1),'\n','\n')
		cat('Relative eigenvalues (% variance in Y; semipartial R-square)','\n')
		cat(EigvalPercent,'\n','\n')	
		cat('Cumulative % variance of species data','\n')
		cat(cumulVar,'\n','\n')
		if(k == 1) {
			cat('NOTE: Since there is a single canonical axis,',
			'the first residual PCA axis','\n')
			cat('has been added to the tables to allow you to',
			'draw a biplot','\n','\n')
			}
		}
	
	# Compute the table of 
	# "Cumulative fit per species as fraction of variance of species" 
	# Canoco manual (1998), p. 174: "With covariables in the analysis, VAR(y) 
	# is unchanged. All fractions are therefore with respect to the original
	# variance [of the species]."
	
	cumfit <- cumulfit.rdaTest(YY, n, p, k, U, Z.sc1, F.sc1)
	# Frac <- FractionBySpecies(YY,U.sc2,n,p,k)
	rownames(cumfit$fit.species) <- var.names
	rownames(cumfit$fit.sites) <- sitenames

	# Create the output list containing the following elements:
#	if(k>1) {
	out <- list(VIF=vif.res, eig.values=eig.values/(n-1), U.sc1=U, U.sc2=U.sc2,
		Z.sc1=Z.sc1, Z.sc2=Z.sc2, F.sc1=F.sc1, F.sc2=F.sc2,
		biplotScores1=posX.sc1, biplotScores2=posX.sc2, 
		fit.species=cumfit$fit.species, fit.sites=cumfit$fit.sites,
		ProbFrda=prob, X.mat=X.mat, Rsquare=Rsquare)
#	} else {
#	out <- list(VIF=vif.res, eig.values=eig.values/(n-1), U.sc1=U, U.sc2=U.sc2,
#		Z.sc1=Z.sc1, Z.sc2=Z.sc2, F.sc1=F.sc1, F.sc2=F.sc2,
#		biplotScores1=posX.sc1, biplotScores2=posX.sc2, ProbFrda=prob, 
#		X.mat=X.mat, Rsquare=Rsquare) }
	class(out) <- "rdaTest"
	out
}

cumulfit.rdaTest <- function(YY, n, p, k, U.sc1, Z.sc1, F.sc1)
# Compute the fractions of the response variables' variances (R2) explained
# by PCA or canonical axes 1, 2, 3, ..., as well as the cumulative fit of 
# the sites along these same axes.
#
# Pierre Legendre, August 2012
{
### Internal function
	sq.length <- function(vec) sum(vec^2)
### End internal function
#
# Extract rownames, etc.
	site.names <- rownames(F.sc1)
	spec.names <- rownames(U.sc1)
#
# Cumulative fit of the species
	Rsq.mat <- matrix(NA,p,k)
	rownames(Rsq.mat) <- spec.names
	colnames(Rsq.mat) <- paste("Axis",1:k,sep="")
	Rsq.mat <- cor(YY, Z.sc1)^2
#
# Cumulative ﬁt of the objects
	cum.obj <- matrix(NA,n,k)
	ref.obj <- apply(YY,1,sq.length)
	cum.obj[,1] <- F.sc1[,1]^2
	if(k>1) for(j in 2:k) cum.obj[,j] <- apply(F.sc1[,1:j],1,sq.length)
	cum.obj <- diag(1/ref.obj) %*% cum.obj
	rownames(cum.obj) <- site.names
	colnames(cum.obj) <- paste("Axis",1:k,sep="")
#
	list(fit.species=t(apply(Rsq.mat,1,cumsum)), fit.sites=cum.obj)
}

probFrda <- function(YY,X,n,p,m,mq,nperm,qr.X,qr.W,qr.XW,SS.Y,SS.Yfit.X,
            SS.Yfit.XW,covar, silent=FALSE, seed, verbose)
# This function carries out a permutation test for the bimultivariate R-square 
# coefficient by permutation of the raw data if there are no covariables, or
# by permutation of residuals of the null model in the presence of covariables.
# Reference: Numerical Ecology, Chapter 11 (Legendre & Legendre 1998).
#
# SS.Y = SS of YY.                                      This is always [a+b+c+d]
# SS.Yfit.X = SS of fitted values of f(Y|X). If covariables present, this is [a]
# SS.Yfit.XW = SS of fitted values of f(Y|(X+W)). If covariables present:[a+b+c]
# 
{
epsilon = .Machine$double.eps
df1 = m
df2 = n-mq-1
if(df2 > 0) {
	if(covar==FALSE) {
		# cat('Debug:  SS.Y =',SS.Y =,'  SS.Yfit.X =',SS.Yfit.X,'\n')
		F.ref = (SS.Yfit.X*df2)/((SS.Y-SS.Yfit.X)*df1)
		} else {
		# cat('Debug: SS.Y =',SS.Y,' SS.Yfit.X =',SS.Yfit.X,' SS.Yfit.XW =', 
		# SS.Yfit.XW,'\n')
		F.ref = (SS.Yfit.X*df2)/((SS.Y-SS.Yfit.XW)*df1)   # [a]/[d]
		Yfit.W = qr.fitted(qr.W, YY)                      # [b + c]
		Yres.W = qr.resid(qr.W, YY)       # [a+b+c+d]-[b+c] = [a+d]
		}
    set.seed(seed)
	nPGE=1
	vec=c(1:n)
	if(covar==FALSE) {
		if(!silent) 
			cat('\n','Test results, permutation of raw data Y','\n',sep="")
		for(i in 1:nperm)
 	  	{
 		# Y.Perm is obtained as follows: YPerm = YY[sample(vec),]
	  	YhatPerm = qr.fitted(qr.X, YY[sample(vec),])
	  	SS.YhatPerm = sum(YhatPerm^2)
 	  	F.perm=(SS.YhatPerm*df2)/((SS.Y-SS.YhatPerm)*df1)
 	  	if(F.perm >= (F.ref-epsilon)) nPGE=nPGE+1
 	  	if(verbose) cat("F =",F.ref,"  F.perm =",F.perm,'\n')
  	  	}
  	  }
  	  else {
  	  # Permute residuals of null model: permute [a+d]. 
		if(!silent) 
			cat('\n','Test results, permutation of residuals of null model',
			'\n',sep="")
		for(i in 1:nperm)
 	  	{
        Yperm = Yfit.W + Yres.W[sample(vec,n),]
        Y.perm = scale(Yperm, center=TRUE, scale=FALSE)
        SS.Yperm = sum(Y.perm^2)               # SS.Yperm is not equal to SS.Y
 		Yperm.hat.W = qr.fitted(qr.W, Yperm)
 		Yperm.hat.XW = qr.fitted(qr.XW, Yperm)
 		SS.Yperm.hat.W = sum(Yperm.hat.W^2)
 		SS.Yperm.hat.XW = sum(Yperm.hat.XW^2)
 	  	F.perm = ((SS.Yperm.hat.XW - SS.Yperm.hat.W) * df2)/((SS.Yperm - 
 	  		SS.Yperm.hat.XW) * df1)
 	  	if(F.perm >= (F.ref-epsilon)) nPGE=nPGE+1
 	  	if(verbose) cat("F =",F.ref,"  F.perm =",F.perm,'\n')
  	  	}
  	  }
	P=nPGE/(nperm+1)
} else {
F.ref = NA
P = NA
if(!silent) cat("\nF-test not computed: residual d.f. <= 0\n")
}

if(verbose) cat('\n')
if(!silent) cat('F =',F.ref,'  Prob(',nperm,'permutations) =',P,'\n')
	
return(list(F=F.ref,nperm=nperm,Prob=P))
}

FractionBySpecies <- function(mat1,U2,n,p,k) 

# This function computes the fractions of the response variables' variances (R2)
# explained by canonical axes 1, 2, 3, ... and by the whole canonical analysis.
#
# mat1 (nxp) is the site-by-species data table Y
# sp.var (p) is a vector containing the species variances
# U2 (pxk) contains the species scores (matrix U.sc2) from PL's rda
# res (pxk) contains the results, found in the output element $rdaFitSpe
# The output element  $VarExpl  contains the % of each species'
#    variance explained by the canonical analysis
#
# Use of the function:
#	Frac=FractionBySpecies(Y.mat,toto$rdaUSc2)
#	Frac$rdaFitSpe   table of Cumulative fit per species as fraction of
#                    variance of species
#	Frac$VarExpl     vector of total % fit per species after all axes
{
	sp.var <- diag(var(mat1))
	res <- matrix(NA,p,k)
	VarExpl <- vector(mode="numeric",p)
	for(i in 1:p) {
#	   var.expl=0.0
#	   for(j in 1:k) {
#	      var.expl=var.expl+(U2[i,j]^2)
#	      res[i,j]=100*var.expl/sp.var[i]
#	   } 
	   res[i,] <- 100*cumsum(U2[i,1:k]^2)/sp.var[i]
	   VarExpl[i] <- res[i,k]
	}
	varnames <- colnames(mat1,do.NULL = FALSE, prefix = "p")
	names(VarExpl) <- varnames
	rownames(res) <- varnames
	colnames(res) <- colnames(res,do.NULL=FALSE,prefix="Cum.axis")
	return(list(rdaFitSpe=res,VarExpl=VarExpl))
}

vif <- function(mat, print.vif=TRUE)
#
# This function returns a vector containing Variance Inflation Factors (VIF)
# for the explanatory variables X. VIFs are the diagonal terms of the inverse
# of the correlation matrix.
#
# Reference: p. 386 of
#    Neter, J. et al. 1996. Applied linear statistical models. 4th edition.
#    Richard D. Irwin Inc., Chicago.
#
# mat = response data matrix to be filtered for collinearity.	
#
# Output:
#
# - The variables that are completely collinear with the previously entered 
#   variables receive a 0 code.
# - The other variables receive Variance Inflation Factors (VIF). 
#   VIF > 10 (or > 20 for other authors) represent high collinearity.
#
# Sebastien Durand and Pierre Legendre, March 2005
{
	threshold <- .Machine$double.eps
	# threshold = value defining the threshold for considering that a covariance 
	# matrix has a zero determinant. A covariance or correlation matrix having a 			
	# zero determinant contains at least one variable that is completely 
	# collinear with the others.
	mat.cor <- cor(mat)
	mm <- FALSE
	for(i in 2:ncol(mat)) {
		look <- c(1:i)
		if(mm) look <- look[-mrk]
		if( det(mat.cor[look,look]) < threshold ) {
			if(print.vif) {
				cat('\n')
				cat("Variable '",colnames(mat)[i],
				"' is collinear: determinant is ", det(mat.cor),'\n')
				}
			if(!mm) {
				mrk <- i
				mm  <- TRUE
			} else {
				mrk <- c(mrk,i)
			}
		}	
	}
	
	# Assign 0 to VIF of variables that are completely collinear with the 
	# previous variables. Compute VIF for the remaining variables
	if(mm) {
		if(ncol(as.matrix(mat[,-mrk])) == 1) {
			vif1 <- 1
			} else {
			vif1 <- diag(solve(cor(mat[,-mrk])))
			}
		vif0 <- rep(0,length(mrk))
		vif  <- c(vif0,vif1)
		nn   <- 1:ncol(mat)
		vif  <- round(vif[sort(c(mrk,nn[-mrk]),index.return=TRUE)$ix],digit=2)
	} else {
		vif  <- round(diag(solve(cor(mat))),digit=2)
	}
	names(vif) <- colnames(mat)
	
	#
	return(vif)
}
