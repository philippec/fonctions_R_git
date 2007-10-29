forward.sel.par <- function(Y, X, alpha = 0.05, K = nrow(X)-1, R2thresh = 0.99, R2more = 0.001, adjR2thresh = 0.99, Yscale = FALSE, verbose=TRUE)
##
## Parametric forward selection of explanatory variables in regression and RDA.
## Y is the response, X is the table of explanatory variables.
##
## If Y is univariate, this function implements FS in regression.
## If Y is multivariate, this function implements FS using the F-test described 
## by Miller and Farr (1971). This test requires that
##   -- the Y variables be standardized,
##   -- the error in the response variables be normally distributed (to be verified by the user).
##
## This function uses 'simpleRDA2' and 'RsquareAdj' developed for 'varpart' in 'vegan'.
##
##                Pierre Legendre & Guillaume Blanchet, May 2007
##
## Arguments --
##
## Y         Response data matrix with n rows and m columns containing quantitative variables.
## X         Explanatory data matrix with n rows and p columns containing quantitative variables.
## alpha     Significance level. Stop the forward selection procedure if the p-value of a variable is higher than alpha. The default is 0.05.
## K         Maximum number of variables to be selected. The default is one minus the number of rows.
## R2thresh  Stop the forward selection procedure if the R-square of the model exceeds the stated value. This parameter can vary from 0.001 to 1.
## R2more    Stop the forward selection procedure if the difference in model R-square with the previous step is lower than R2more. The default setting is 0.001.
## adjR2thresh Stop the forward selection procedure if the adjusted R-square of the model exceeds the stated value. This parameter can take any value (positive or negative) smaller than 1.
## Yscale    Standardize the variables in table Y to variance 1. The default setting is FALSE. The setting is automatically changed to TRUE if Y contains more than one variable. This is a validity condition for the parametric test of significance (Miller and Farr 1971).
##
## Reference:
## Miller, J. K., and S. D. Farr. 1971. Bimultivariate redundancy: a comprehensive measure of 
##    interbattery relationship. Multivariate Behavioral Research 6: 313-324.

{
  require(vegan)
  FPval <- function(R2cum,R2prev,n,mm,p)
    ## Compute the partial F and p-value after adding a single explanatory variable to the model.
    ## In FS, the number of df of the numerator of F is always 1. See Sokal & Rohlf 1995, eq 16.14.
    ## 
    ## The amendment, based on Miller and Farr (1971), consists in multiplying the numerator and  
    ## denominator df by 'p', the number of variables in Y, when computing the p-value.
    ##
    ##                Pierre Legendre, May 2007
    {
      df2 <- (n-1-mm)
      Fstat <- ((R2cum-R2prev)*df2) / (1-R2cum)
      pval <- pf(Fstat,1*p,df2*p,lower.tail=FALSE)
      return(list(Fstat=Fstat,pval=pval))
    }

  Y <- as.matrix(Y)
  X <- apply(as.matrix(X),2,scale,center=TRUE,scale=TRUE)
  var.names = colnames(as.data.frame(X))
  n <- nrow(X)
  m <- ncol(X)
  if(nrow(Y) != n) stop("Numbers of rows not the same in Y and X")
  p <- ncol(Y)
  if(p > 1) {
    Yscale = TRUE
    if(verbose) cat("The variables in response matrix Y have been standardized",'\n')
  }
  Y <- apply(Y,2,scale,center=TRUE,scale=Yscale)
  SS.Y <- sum(Y^2)
  
  X.out <- c(1:m)
  
  ## Find the first variable X to include in the model
  R2prev <- 0
  R2cum <- 0
  for(j in 1:m) {
    toto <- simpleRDA2(Y,X[,j],SS.Y)
    if(toto$Rsquare > R2cum) {
      R2cum <- toto$Rsquare
      no.sup <- j
    }
  }
  mm <- 1
  FP <- FPval(R2cum,R2prev,n,mm,p)
  if(FP$pval <= alpha) {
    adjRsq <- RsquareAdj(R2cum,n,mm)
    res1 <- var.names[no.sup]
    res2 <- no.sup
    res3 <- R2cum
    res4 <- R2cum
    res5 <- adjRsq
    res6 <- FP$Fstat
    res7 <- FP$pval
    X.out[no.sup] <- 0
    delta <- R2cum
  } else {
    stop("Procedure stopped (alpha criterion): pvalue for variable ",no.sup," is ",FP$pval)
  }
  
  ## Add variables X to the model
  while((FP$pval <= alpha) & (mm <= K) & (R2cum <= R2thresh) & (delta >= R2more) & (adjRsq <= adjR2thresh)) {
    mm <- mm+1
    R2prev <- R2cum
    R2cum <- 0
    for(j in 1:m) {
      if(X.out[j] != 0) {
        toto <- simpleRDA2(Y,X[,c(res2,j)],SS.Y)
        if(toto$Rsquare > R2cum) {
          R2cum <- toto$Rsquare
          no.sup <- j
        }
      }
    }
    FP <- FPval(R2cum,R2prev,n,mm,p)
    delta <- R2cum-R2prev
    adjRsq <- RsquareAdj(R2cum,n,mm)
    res1 <- c(res1,var.names[no.sup])
    res2 <- c(res2,no.sup)
    res3 <- c(res3,delta)
    res4 <- c(res4,R2cum)
    res5 <- c(res5,adjRsq)
    res6 <- c(res6,FP$Fstat)
    res7 <- c(res7,FP$pval)
    X.out[no.sup] <- 0
  }
  if(verbose) {
    if(FP$pval > alpha)  cat("Procedure stopped (alpha criterion): pvalue for variable ",no.sup," is ",FP$pval,'\n') 
    if(mm > K)           cat("Procedure stopped (K criterion): mm = ",mm," is larger than ",K," after including variable ",no.sup,'\n') 
    if(R2cum > R2thresh) cat("Procedure stopped (R2thresh criterion): R2cum for variable ",no.sup," is ",R2cum,'\n') 
    if(delta < R2more)   cat("Procedure stopped (R2more criterion): delta for variable ",no.sup," is ",delta,'\n') 
    if(adjRsq>adjR2thresh) cat("Procedure stopped (adjR2thresh criterion): adjRsq for variable ",no.sup," is ",adjRsq,'\n') 
  }
  
  res <- data.frame(res1,res2,res3,res4,res5,res6,res7)
  colnames(res) <- c("variable","order","R2","R2cum","AdjR2Cum","F","pval")
  if((FP$pval > alpha) | (mm > K) | (R2cum > R2thresh) | (delta < R2more) | (adjRsq > adjR2thresh))  res <- res[1:(mm-1),]
  
  return(res)
}

