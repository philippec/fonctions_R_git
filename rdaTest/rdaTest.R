rdaTest <- function(YY.mat, XX.mat, WW.mat=NULL, 
           scale.Y=FALSE, testF=NULL, nperm=NULL, print.results=TRUE, print.cum=FALSE)
#
# Function RDA for simple or partial RDA, with permutation tests,
# computed as described in Numerical Ecology, Chapter 11 (Legendre & Legendre 1998)
#
# Pierre Legendre, UniversitŽ de MontrŽal, March-April 2005, May-August 2007
#
# Modified June and August 2006: This version runs correctly with a single   
# explanatory variable in XX and a single response variable in YY.
# Since there is a single canonical axis in either case, the first PCA axis of
# the residuals is added to the tables to allow the drawing of biplots.
#
#
# PARAMETERS:
#
# YY.mat (nxp) is the site-by-species data table
# XX.mat (nxm) contains the explanatory variables. The number of variables
#    'm' is computed after eliminating collinear explanatory variables, if any
# WW.mat (nxm) contains the covariables, if any
#
# scale.Y contains a logical value: should YY.mat be standardized, or not
# testF: when NULL, the program will ask the user if he/she wishes to test the  
#    F statistic. If testF is TRUE or FALSE, no question will be asked; 
#    the program will perform the test, or not, in accordance with that indication.
# nperm: number of permutation for the F test. If NULL, a question will be asked.
# print.results: prints the main rdaTest results on the screen.
# print.cum: prints the fractions of the response variable's (e.g. species) variances
#    explained by canonical axes 1, 2, 3, ... and by the whole canonical analysis.
#
#
# The function first returns a list of XX.mat variables with null variances, if any.
# The function then returns an output list containing the following ELEMENTS:
#
# return(list(VIF=vif.res, canEigval=canEigval, U=U, USc2=UL, F=F, Z=Z, FSc2=FSc2, 
#        ZSc2=ZSc2, biplotScores1=posX, biplotScores2=posXSc2, FitSpe=Frac$rdaFitSpe, 
#        VarExpl=Frac$VarExpl, X.mat=X.mat))
#
# VIF: variance inflation factors for the explanatory variables X; 
#    the value is 0 for entirely collinear variables
# canEigval: canonical eigenvalues
# U (pxk): canonical eigenvectors normalized to 1 (scaling 1)
# USc2 (pxk): canonical eigenvectors normalized to sqrt(eigenvalue) (scaling 2)
# F (nxk): matrix of object scores (scaling 1)
# Z (nxk): matrix of fitted object scores (scaling 1)
# FSc2: matrix of object scores (scaling 2)
# ZSc2: matrix of fitted object scores (scaling 2)
# biplotScores1: biplot scores of explanatory variables (scaling 1)
# biplotScores2: biplot scores of explanatory variables (scaling 2)
# FitSpe: table of cumulative fit per species (in %) as fraction of variance of species
# VarExpl: vector of total % fit per species after all canonical axes 
# X.mat: original X matrix (required by the plotting function)
#
#
# Use of the function 'rdaTest':
#
# First set of examples: use all default values, except the data file names
# result <- rdaTest(Y, X)         # example of simple RDA
# result <- rdaTest(Y, X, W)      # example of partial RDA
#
# Second set of examples: ask for the permutation test in the function parameters
# result <- rdaTest(Y, X, testF=TRUE, nperm=999)
# result <- rdaTest(Y, X, W, testF=TRUE, nperm=999)
#
# How to obtain the non-canonical axes:
# result.noncan <- rdaTest(Y, Y, X)           # Example with explanatory variables X
# result.noncan <- rdaTest(Y, Y, cbind(X,W))  # Example with X and covariables W
#
# 345678901234567890123456789012345678901234567890123456789012345678901234567890
{
	library(MASS)
	if(is.logical(scale.Y)){
	}else{
		stop("Wrong operator; 'center.Y' should be either 'FALSE' or 'TRUE'")	
	}
		
	# Read the data tables. Transform them into matrices Y, X, and W
	Y.mat=as.matrix(YY.mat)
	X.mat=as.matrix(XX.mat)
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
		vif.res=vif(X,print.vif=print.results)
		if(length(which(vif.res==0))!=0) X.cor = X[,-which(vif.res==0)]
		} else { vif.res=rep(1,m1) }
	m=ncol(X.cor)
	
	if(print.results==TRUE & m > 1){ 
		cat('\n')
		cat('Variance Inflation Factors (VIF)','\n')
		print(vif.res)
	}		

	qq=0
	mq=m
	# Is there a matrix W containing covariables?
	if(length(WW.mat)!=0){
		covar=TRUE
		W.mat=as.matrix(WW.mat)
		qq=ncol(W.mat)
		W=apply(W.mat,2,scale,center=TRUE,scale=TRUE)
    	if(print.results==TRUE) cat("\nThere is a covariance matrix\n")
    	# Find the rank of W using QR decomposition
    	QR.W = qr(W, tol = 1e-06)
    	q = QR.W$rank
        # Find the rank of cbind(X,W) using QR decomposition
    	QR.XW = qr(cbind(X,W), tol = 1e-06)
    	mq = QR.XW$rank
	}
	else{
		covar=FALSE
    	if(print.results==TRUE) cat("\nThere is no covariance matrix\n")
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
		adjRsq = 1-((1-Rsquare)*totalDF/residualDF)
	}

	# Test significance of the canonical F statistic
	testFF="N"
	if(is.null(testF)){
		cat('\n')
		cat("Test the bimultivariate F-statistic? Type 'Y' if yes, 'N' if no.",'\n')
		testFF=toupper(scan(file="",what="character",nlines=1,quiet=T))
	    if(testFF=="y") testFF="Y"
	} else { if(testF==TRUE) testFF="Y" }

	if(testFF=="Y") {
		if(is.null(nperm)){
			cat('\n','How many permutations? Ex. 499, 999, ...','\n',sep="")
			nper=toupper(scan(file="",what="integer",nlines=1,quiet=T))
			nperm=as.integer(nper)
		}
	}
	if(print.results==TRUE){
   		cat('\n','----------','\n','\n',sep="")
		cat('Bimultivariate redundancy statistic (canonical R-square):','\n')
		cat('R-square =',Rsquare)
		if(covar==FALSE) cat(';   adjusted R-square = ',adjRsq)
		cat('\n')
	}

	if(testFF=="Y") {prob<-probFrda(Y,X,n,p,m,mq,nperm,projX,projW,projXW,SS.Y,SS.Yfit.X,
	                      covar,SS.Yfit.XW,print.results=print.results)
	}else{ prob = "Not tested" }
	
	# PCA portion of RDA: eigenanalysis, then compute the F and Z matrices, etc.
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
	F = Y %*% U
	Z = Yfit.X %*% U
	
	# for (scaling==2)
	lambdaSc2=vector(mode="numeric",k)
	for(j in 1:k) {
   		lambdaSc2[j]=sqrt(Yhat.val[j])
   	}
	if(k == 1) {
		UL=U * lambdaSc2[1]
   		FSc2 = F * (1/lambdaSc2[1])
   		ZSc2 = Z * (1/lambdaSc2[1])
   	}
	if(k > 1) {
		UL=U %*% diag(lambdaSc2)
		FSc2 = F %*% diag(1/lambdaSc2)
		ZSc2 = Z %*% diag(1/lambdaSc2)
   	}
	# Compute the 'Biplot scores of environmental variables' --
	# First, compute the correlations between X and Z, or X and ZSc2
	#
	# Positions for scaling 1: compute the correlations...
	corXZ=cor(X,Z)
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
	posXSc2=corXZ   # which is cor(X,Z); posXSc2=cor(X,ZSc2) produces the same result
	
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
       F = cbind(F,Yres.F[,1])
       Z = cbind(Z,Yres.F[,1])
       UL = cbind(UL,Yres.U2[,1])
       FSc2 = cbind(FSc2,Yres.G[,1])
       ZSc2 = cbind(ZSc2,Yres.G[,1])
       posX = cbind(posX,0)
       posXSc2 = cbind(posXSc2,0)
	   }
	
	# Print results
   	if(print.results==TRUE){
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
	} else {}
	
	# Compute the "Cumulative fit per species as fraction of variance of species" table
	Frac=FractionBySpecies(Y.mat,UL,n,p,k,print.cum=print.cum)

	# Create the output list containing the following elements:
	#
	# VIF=vif.res: Variance Inflation Factors (VIF)
	# canEigval=canEigval: vector of canonical eigenvalues
	# U=U: canonical eigenvectors normalized to 1 (scaling 1)
	# USc2=UL: canonical eigenvectors normalized to sqrt(eigenvalue) (scaling 2)
	# F=F: matrix of object scores (scaling 1)
	# Z=Z: matrix of fitted object scores (scaling 1)
	# FSc2=FSc2: matrix of object scores (scaling 2)
	# ZSc2=ZSc2: matrix of fitted object scores (scaling 2)
	# biplotScores1=posX: biplot scores of explanatory variables (scaling 1)
	# biplotScores2=posXSc2: biplot scores of explanatory variables (scaling 2)
	# FitSpe=Frac$rdaFitSpe: table of cumulative fit per species (in %) 
	#    as fraction of variance of species
	# VarExpl=Frac$VarExpl: vector of total % fit per species after all canonical axes 
	# ProbFrda=prob: probability associated with F test of the canonical relationship
	# X.mat=X.mat: original X matrix (required by the plotting function)

	return(list(VIF=vif.res, canEigval=canEigval, U=U, USc2=UL, F=F, Z=Z, FSc2=FSc2, 
	       ZSc2=ZSc2, biplotScores1=posX, biplotScores2=posXSc2, FitSpe=Frac$rdaFitSpe, 
	       VarExpl=Frac$VarExpl, ProbFrda=prob, X.mat=X.mat, Rsq=Rsquare))
}


probFrda <- function(Y,X,n,p,m,mq,nperm,projX,projW,projXW,SS.Y,SS.Yfit.X,covar,
                     SS.Yfit.XW,print.results=TRUE)
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
    epsilon=1e-15
	df1=m
	df2=n-mq-1
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
		if(print.results==TRUE) 
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
		if(print.results==TRUE) 
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
	if(print.results==TRUE)
		cat('F =',Fref,'  Prob(',nperm,'permutations) =',P,'\n')
	return(list(F=Fref,nperm=nperm,Prob=P))
}


FractionBySpecies <- function(mat1,mat3,n,p,k,print.cum=FALSE) 

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
	if(print.cum==TRUE){
	#	cat('Debug: Values of n =',n,'  p =',p,'  k =',k,'\n','\n')
		cat('Cumulative fit of species as fraction of their variance ',
			'after 1, 2, 3 ... axes','\n')
		print(mat4)
		cat('\n')
	}else{}
	return(list(rdaFitSpe=mat4,VarExpl=VarExpl))
}


vif<-function(mat, threshold=1.0e-12, print.vif=TRUE){
#
# This function returns a vector containing Variance Inflation Factors (VIF)
# for the explanatory variables X. VIFs are the diagonal terms of the inverse
# of the correlation matrix.
#
# Reference: p. 386 of
#    Neter, J. et al. 1996. Applied linear statistical models. 4th edition.
#    Richard D. Irwin Inc., Chicago.
#
# mat = matrix to be filtered for collinearity.	
# threshold = value defining the threshold for considering that a covariance matrix
#    has a zero determinant. A covariance or correlation matrix having a zero determinant
#    contains at least one variable that is completely collinear with the others.
#
# Output:
#
# - The variables that are completely collinear with the previously entered variables receive 0.
# - The other variables receive Variance Inflation Factors (VIF). VIF > 10 (or > 20 for other
#   authors) represent high collinearity.
#
# Sebastien Durand and Pierre Legendre, March 2005
#
	library(MASS)

	# Assign 0 to the variables that are completely collinear with the previous variables
	mat.cor <- cor(mat)
	mm = 0
	for(i in 2:ncol(mat)){
		look<-c(1:i)
		if(mm!=0)
			look <- look[-mrk]
		if( det(mat.cor[look,look]) < threshold ) {
			if(print.vif==T) {
				cat('\n')
				cat("Variable '",colnames(mat)[i],"' is collinear: determinant is ", 
				    det(mat.cor),'\n')
			}
			if(mm==0) {
				mrk <- i
				mm <- 1
			} else {
				mrk <- c(mrk,i)
			}
		}	
	}
	
	# Compute VIF for the remaining variables
	if(mm==1){
		vif1<-diag(ginv(cor(mat[,-mrk])))
		vif0<-rep(0,length(mrk))
		vif<-c(vif0,vif1)
		nn<-1:ncol(mat)
		vif <- round(vif[sort(c(mrk,nn[-mrk]),index.return=T)$ix],digit=2)
	}else{
		vif <- round(diag(ginv(cor(mat))),digit=2)
	}
	names(vif)<-colnames(mat)
	
	#
	return(vif)
}
