'rdaTest' <- 
	function(Y, X, W=NULL, scale.Y=FALSE, test.F=NULL, nperm=NULL, silent=FALSE)
{
	if(!is.logical(scale.Y)) stop("Wrong operator; 'center.Y' should be either 'FALSE' or 'TRUE'")	
		
	# Read the data tables. Transform them into matrices Y, X, and W
	Y.mat=as.matrix(Y)
	X.mat=as.matrix(X)
	sitenames<-rownames(X.mat)
	if(is.null(rownames(X.mat)) & is.null(rownames(Y.mat))){
		sitenames<-paste("Site",1:nrow(X.mat),sep="")
	}else{
		if(is.null(rownames(X.mat))){
			sitenames<-rownames(Y.mat)
		}else{
			sitenames<-rownames(X.mat)
		}
	}
	
	# If(scale.Y == TRUE), check the presence of Y.mat variables with null variances
	if(scale.Y == TRUE) {
	mat.sd <- sd(Y.mat)
	problem = FALSE
	for(i in 1:length(mat.sd)) {
	   if(mat.sd[i] == 0) {
	      cat(" Variable No.",i," in Y has a null variance",'\n')
	      problem = TRUE
	      }
	   }
	if(problem == TRUE) stop("The program was stopped. Verify/modify your Y matrix.")	
	}

	# Check the presence of X.mat variables with null variances
	mat.sd <- sd(X.mat)
	problem = FALSE
	for(i in 1:length(mat.sd)) {
	   if(mat.sd[i] == 0) {
	      cat(" Variable No.",i," in X has a null variance",'\n')
	      problem = TRUE
	      }
	   }
	if(problem == TRUE) stop("The program was stopped. Verify/modify your X matrix.")
	
	# Pre-treatment: centre or standardize Y, standardize X
	n=nrow(Y.mat)
	p=ncol(Y.mat)
	Y = apply(Y.mat,2,scale,center=TRUE,scale=scale.Y)
	rownames(Y)<-sitenames
	X = apply(X.mat,2,scale,center=TRUE,scale=TRUE)
	
	# Identify collinear columns in X and remove them
	m1 = ncol(X)
	X.cor = X
	if(m1 > 1) {
		vif.res = vif(X,print.vif=!silent)
		if(length(which(vif.res==0))!=0) X.cor = X[,-which(vif.res==0)]
		} else { vif.res=rep(1,m1) }
	m=ncol(as.matrix(X.cor))
	
	if((!silent) & m > 1){ 
		cat('\n')
		cat('Variance Inflation Factors (VIF)','\n')
		print(vif.res)
	}		

	qq=0
	mq=m
	# Is there a matrix W containing covariables?
	if(length(W)!=0){
		covar=TRUE
		W.mat=as.matrix(W)
		qq=ncol(W.mat)
		W=apply(W.mat,2,scale,center=TRUE,scale=TRUE)
    	if(!silent) cat("\nThere is a matrix of covariables\n")
    	# Find the rank of W using QR decomposition
    	QR.W = qr(W, tol = 1e-06)
    	q = QR.W$rank
        # Find the rank of cbind(X,W) using QR decomposition
    	QR.XW = qr(cbind(X,W), tol = 1e-06)
    	mq = QR.XW$rank
	}
	else{
		covar=FALSE
    	if(!silent) cat("\nThere is no matrix of covariables\n")
	}

	# If covariables W are present, regress X on W. Obtain X.res
	if(covar==TRUE) {
		invW = ginv(t(W) %*% W)
		projW = W %*% invW %*% t(W)
		X.res = X.cor - projW %*% X.cor
		
		XW=cbind(X.cor,W)
		invXW = ginv(t(XW) %*% XW)
		projXW = XW %*% invXW %*% t(XW)
		XW.fit = projXW %*% Y
		SS.Yfit.XW = sum(XW.fit*XW.fit)
	}else{
		X.res = X.cor
	}
	
	# Compute projector of X and Yfit.X
	invX = ginv(t(X.res) %*% X.res)
	projX = X.res %*% invX %*% t(X.res)
	Yfit.X = projX %*% Y
	rownames(Yfit.X)<-sitenames
	
	# Compute R-square and adjusted R-square
	SS.Y = sum(Y^2)
	SS.Yfit.X = sum(Yfit.X^2)
	
	# SS.Y = sum(diag(cov(Y)))
	# SS.Yfit.X = sum(diag(cov(Yfit.X)))
	Rsquare = SS.Yfit.X/SS.Y
	if(covar==FALSE) {
		totalDF=n-1
		residualDF=n-m-1
		if(residualDF > 0) {
			adjRsq = 1-((1-Rsquare)*totalDF/residualDF)
			} else {
			adjRsq = NA
			if(!silent) cat("\nAdjusted R-square not computed: residual d.f. <= 0\n")
			}
	}

	# Test significance of the canonical F statistic
	testFF="N"
	if(is.null(test.F)){
		cat('\n')
		cat("Test the bimultivariate F-statistic? Type 'Y' if yes, 'N' if no.",'\n')
		testFF=toupper(scan(file="",what="character",nlines=1,quiet=TRUE))
	    if(testFF=="y") testFF="Y"
	} else { if(test.F==TRUE) testFF="Y" }

	if(testFF=="Y") {
		if(is.null(nperm)){
			cat('\n','How many permutations? Ex. 499, 999, ...','\n',sep="")
			nper=toupper(scan(file="",what="integer",nlines=1,quiet=TRUE))
			nperm=as.integer(nper)
		}
	}
	if(!silent){
   		cat('\n','----------','\n','\n',sep="")
		cat('Bimultivariate redundancy statistic (canonical R-square):','\n')
		cat('R-square =',Rsquare)
		if(covar==FALSE) cat(';   adjusted R-square = ',adjRsq)
		cat('\n')
	}

	if(testFF=="Y") {
	prob<-probFrda(Y,X,n,p,m,mq,nperm,projX,projW,projXW,SS.Y,SS.Yfit.X,
	                      covar,SS.Yfit.XW,silent=silent)
	} else { prob = "Not tested" }
	
	# PCA portion of RDA: eigenanalysis, then compute the F.mat and Z.mat matrices, etc.
	SS.Y=SS.Y/(n-1)
	Yhat.cov = cov(Yfit.X)
	Yhat.eig = eigen(Yhat.cov)
	Yhat.val=Yhat.eig$values
	axenames<-paste("Axis",1:length(Yhat.val),sep="")
	names(Yhat.val)<-axenames
	Yhat.vec=Yhat.eig$vectors
	rownom<-colnames(X.res,do.NULL = FALSE, prefix = "m")
	colnames(Yhat.vec)<-axenames
	rownames(Yhat.vec)<-colnames(Y,do.NULL = FALSE, prefix = "p")
	# How many canonical eigenvalues?
	kk=min((n-1),p,m)
	k = length(which(Yhat.val > 0.00000001))
	# k=0
	# for(i in 1:kk) { if(Yhat.val[i] > 0.00000001) k=k+1 }
	#
	canEigval=Yhat.val[1:k]
	EigvalCanoco = 100*canEigval/SS.Y
	cumulVar=vector(mode="numeric",k)
	cumulVar[1]=EigvalCanoco[1]
	if(k > 1) {
		for(j in 2:k) { cumulVar[j]=cumulVar[j-1]+EigvalCanoco[j] } 
   	}
	
	# for (scaling==1)
	U = as.matrix(Yhat.vec[,1:k])
	F.mat = Y %*% U
	Z.mat = Yfit.X %*% U
	
	# for (scaling==2)
	lambdaSc2=vector(mode="numeric",k)
	for(j in 1:k) {
   		lambdaSc2[j]=sqrt(Yhat.val[j])
   	}
	if(k == 1) {
		UL=U * lambdaSc2[1]
   		FSc2 = F.mat * (1/lambdaSc2[1])
   		ZSc2 = Z.mat * (1/lambdaSc2[1])
   	}
	if(k > 1) {
		UL=U %*% diag(lambdaSc2)
		FSc2 = F.mat %*% diag(1/lambdaSc2)
		ZSc2 = Z.mat %*% diag(1/lambdaSc2)
   	}
	# Compute the 'Biplot scores of environmental variables' --
	# First, compute the correlations between X and Z.mat, or X and ZSc2
	#
	# Positions for scaling 1: compute the correlations...
	corXZ=cor(X,Z.mat)
	# ...then, weigh these correlations
	# by the diagonal matrix 'D' of weights  sqrt(lambda(k)/SS.Y)
	if(k == 1) { 
	   D = as.matrix(sqrt(Yhat.eig$values[1]/SS.Y))
	   } else {
   	   D = diag(sqrt(Yhat.eig$values[1:k]/SS.Y))
   	   }
   	posX = corXZ %*% D
   	#
	# Positions for scaling 2:
	posXSc2=corXZ   # which is cor(X,Z.mat); posXSc2=cor(X,ZSc2) produces the same result
	
	# If k = 1, compute the first PCA axis of the residuals and add it to all tables
	if(k == 1) {
	   Yres.X = Y - Yfit.X
       Yres.cov = cov(Yres.X)
       Yres.eig = eigen(Yres.cov)
       Yres.U = Yres.eig$vectors
       Yres.F = Yres.X %*% Yres.U
	   if(p == 1) {
          Yres.U2 = Yres.U %*% Yres.eig$value^(0.5)
          Yres.G = Yres.F %*% Yres.eig$value^(-0.5)
   	      } else {
          Yres.U2 = Yres.U %*% diag(Yres.eig$value^(0.5))
          Yres.G = Yres.F %*% diag(Yres.eig$value^(-0.5))
          }

       canEigval = c(canEigval,Yres.eig$values[1])
       U = cbind(U,Yres.U[,1])
       F.mat = cbind(F.mat,Yres.F[,1])
       Z.mat = cbind(Z.mat,Yres.F[,1])
       UL = cbind(UL,Yres.U2[,1])
       FSc2 = cbind(FSc2,Yres.G[,1])
       ZSc2 = cbind(ZSc2,Yres.G[,1])
       posX = cbind(posX,0)
       posXSc2 = cbind(posXSc2,0)
	   }
	
	# Print results
   	if(!silent){
   		cat('\n','Number of objects:  n = ',n,'\n',sep="")
   		cat('Number of response variables in Y:  p =',p,'\n')
   		cat('Number of explanatory variables in X: ',m1,';  Rank of X:  m =',m,'\n')
   		if(covar==TRUE) {
   		   cat('Number of covariables in W: ',qq,';  Rank of W:  q =',q,'\n')
   		   cat('Rank of cbind(X,W):  mq =',mq,'\n')
   		   }
   		cat('Number of canonical eigenvalues:  k =',k,'\n')
   		cat('\n')
		cat('Total variance =',SS.Y,'\n','\n')
		cat('Eigenvalues','\n')
		cat(canEigval,'\n','\n')
		cat('Relative eigenvalues (% variance)','\n')
		cat(EigvalCanoco,'\n','\n')	
		cat('Cumulative % variance of species data','\n')
		cat(cumulVar,'\n','\n')
		if(k == 1) {
		cat('NOTE: Since there is a single canonical axis, the first residual PCA axis','\n')
		cat('has been added to the tables to allow you to draw a biplot','\n','\n')
		}
	}
	
	# Compute the "Cumulative fit per species as fraction of variance of species" table
	Frac=FractionBySpecies(Y,UL,n,p,k)

	# Create the output list containing the following elements:
	out <- list(VIF=vif.res, canEigval=canEigval, U=U, USc2=UL, F.mat=F.mat, 
		Z.mat=Z.mat, FSc2=FSc2, ZSc2=ZSc2, biplotScores1=posX, 
		biplotScores2=posXSc2, FitSpe=Frac$rdaFitSpe, 
		VarExpl=Frac$VarExpl, ProbFrda=prob, X.mat=X.mat, Rsq=Rsquare)
	class(out) <- "rdaTest"
	out
}


probFrda <- function(Y,X,n,p,m,mq,nperm,projX,projW,projXW,SS.Y,SS.Yfit.X,covar,
                     SS.Yfit.XW,silent=FALSE)
# This function carries out a permutation test for the bimultivariate R-square 
# coefficient by permutation of the raw data if there are no covariables, or
# by permutation of the residuals of the null model in the presence of covariables.
# Reference: Numerical Ecology, Chapter 11 (Legendre & Legendre 1998).
#
# SS.Y = SS of Y.                            If covariables present, this is [a+b+c+d]
# SS.Yfit.X = SS of fitted values of f(Y|X).       If covariables present, this is [a]
# SS.Yfit.XW = SS of fitted values of f(Y|(X+W)).      If covariables present: [a+b+c]
# 
{
epsilon = .Machine$double.eps
df1=m
df2=n-mq-1
if(df2 > 0) {
	if(covar==FALSE) {
	#	cat('Debug:  SS.Y =',SS.Y,'  SS.Yfit.X =',SS.Yfit.X,'\n')
		Fref=(SS.Yfit.X*df2)/((SS.Y-SS.Yfit.X)*df1)
	}
	else {
	# cat('Debug: SS.Y =',SS.Y,' SS.Yfit.X =',SS.Yfit.X,' SS.Yfit.XW =',SS.Yfit.XW,'\n')
		Fref=(SS.Yfit.X*df2)/((SS.Y-SS.Yfit.XW)*df1)
		Yfit.W = projW %*% Y    # [b+c]
		Yres.W = Y - Yfit.W     # [a+b+c+d] - [b+c] = [a+d]
	}
	nPGE=1
	vec=c(1:n)
	if(covar==FALSE) {
		if(!silent) 
			cat('\n','Test results, permutation of raw data Y','\n',sep="")
		for(i in 1:nperm)
 	  	{
 		# YPerm is obtained as follows: YPerm = Y[sample(vec,n),]
	  	YhatPerm = projX %*% Y[sample(vec,n),]
	  	SS.YhatPerm = sum(YhatPerm*YhatPerm)
 	  	Fper=(SS.YhatPerm*df2)/((SS.Y-SS.YhatPerm)*df1)
 	  	if(Fper >= (Fref-epsilon)) nPGE=nPGE+1
  	  	}
  	  }
  	  else {
  	  # Permute residuals of null model: permute [a+d].  SS.Yperm is not equal to SS.Y
		if(!silent) 
			cat('\n','Test results, permutation of residuals of null model','\n',sep="")
		for(i in 1:nperm)
 	  	{
        Yperm = Yfit.W + Yres.W[sample(vec,n),]
        Y.perm = apply(Yperm,2,scale,center=TRUE,scale=FALSE)
        SS.Yperm = sum(Y.perm*Y.perm)
 		Yperm.hat.W = projW %*% Yperm
 		Yperm.hat.XW = projXW %*% Yperm
 		SS.Yperm.hat.W = sum(Yperm.hat.W*Yperm.hat.W)
 		SS.Yperm.hat.XW = sum(Yperm.hat.XW*Yperm.hat.XW)
 	  	Fper=((SS.Yperm.hat.XW - SS.Yperm.hat.W) * df2)/((SS.Yperm - SS.Yperm.hat.XW) * df1)
 	  	if(Fper >= (Fref-epsilon)) nPGE=nPGE+1
  	  	}
  	  }
	P=nPGE/(nperm+1)
} else {
Fref = NA
P = NA
if(!silent) cat("\nF-test not computed: residual d.f. <= 0\n")
}

if(!silent) cat('F =',Fref,'  Prob(',nperm,'permutations) =',P,'\n')
	
return(list(F=Fref,nperm=nperm,Prob=P))
}


FractionBySpecies <- function(mat1,mat3,n,p,k) 

# This function computes the fraction of the response variable variance
# explained by canonical axes 1, 2, 3, ... and by the whole canonical analysis.
#
# mat1 (nxp) is the site-by-species data table Y
# spSS (p) is a vector containing the species sums-of-squares
# mat3 (pxk) contains the species scores (matrix UL) from PL's rda
# mat4 (pxk) contains the results, found in the output element $rdaFitSpe
# The output element  $VarExpl  contains the % of each species'
#    variance explained by the canonical analysis
#
# Use of the function:
#	Frac=FractionBySpecies(Y.mat,toto$rdaUSc2)
#	Frac$rdaFitSpe   table of Cumulative fit per species as fraction of
#                    variance of species
#	Frac$VarExpl     vector of total % fit per species after all axes
{
	#	n=nrow(mat1)
	#	p=ncol(mat1)
	#	k=ncol(mat3)
	spSS=diag(var(mat1))
	mat4=matrix(NA,p,k)
	VarExpl=vector(mode="numeric",p)
	for(i in 1:p) {
	   ss=0.0
	   for(j in 1:k) {
	      ss=ss+(mat3[i,j]^2)
	      mat4[i,j]=100*ss/spSS[i]
	   } 
	   VarExpl[i]=mat4[i,k]
	}
	varnames<-colnames(mat1,do.NULL = FALSE, prefix = "p")
	names(VarExpl)<-varnames
	rownames(mat4)<-varnames
	colnames(mat4)<-colnames(mat4,do.NULL=FALSE,prefix="Cum.axis")
	return(list(rdaFitSpe=mat4,VarExpl=VarExpl))
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
	library(MASS)
	threshold <- .Machine$double.eps
	# threshold = value defining the threshold for considering that a covariance 
	# matrix has a zero determinant. A covariance or correlation matrix having a 	# zero determinant contains at least one variable that is completely 
	# collinear with the others.
	mat.cor <- cor(mat)
	mm = 0
	for(i in 2:ncol(mat)) {
		look <- c(1:i)
		if(mm != 0)
			look <- look[-mrk]
		if( det(mat.cor[look,look]) < threshold ) {
			if(print.vif) {
				cat('\n')
				cat("Variable '",colnames(mat)[i],"' is collinear: determinant is ", det(mat.cor),'\n')
				}
			if(mm == 0) {
				mrk <- i
				mm  <- 1
			} else {
				mrk <- c(mrk,i)
			}
		}	
	}
	
	# Assign 0 to VIF of variables that are completely collinear with the 
	# previous variables. Compute VIF for the remaining variables
	if(mm==1) {
		if(ncol(as.matrix(mat[,-mrk])) == 1) {
			vif1 <- 1
			} else {
			vif1 <- diag(ginv(cor(mat[,-mrk])))
			}
		vif0 <- rep(0,length(mrk))
		vif  <- c(vif0,vif1)
		nn   <- 1:ncol(mat)
		vif  <- round(vif[sort(c(mrk,nn[-mrk]),index.return=TRUE)$ix],digit=2)
	} else {
		vif  <- round(diag(ginv(cor(mat))),digit=2)
	}
	names(vif) <- colnames(mat)
	
	#
	return(vif)
}
