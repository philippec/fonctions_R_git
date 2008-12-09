quickPCNM <- function(Y,space,truncd=0,method="fwd",PCNM=NULL,alpha=0.05,rangexy=FALSE,detrend=TRUE){

#              *** Quick exploratory PCNM analysis ***
#                  Version 7.7.1 (27 november 2008)
#
#   PCNM analysis for quick assessment of spatial structures, either on
#   the basis of geographical coordinates (which can be one- or two-
#   dimensional) or wich a set of pre-computed complex spatial variables
#   (PCNM, MEM, polynomials...). 
#   If ONLY geographical coordinates are provided, the function automatically
#   computes PCNM variables.
#   Geographical coordinates must be provided even if complex spatial variables
#   are given as well.
#   If PCNM or other complex variables are provided in the PCNM argument, 
#   the function runs the spatial analysis directly on these variables.
#   Be careful NOT to provide enough spatial variables to saturate the model,
#   i.e., the number of spatial variables must be smaller than n-1. Otherwise
#   the global test will be nonsignificant and the analysis will be stopped.
#
#   The RDA is computed on a covariance matrix. If one wishes an RDA on
#   a correlation matrix instead, the dependent variables have to be
#   standardized to zero mean and unit variance prior to the analysis.
#
#   This function offers 4 different methods of selection of explanatory
#   variables (see below).
#
#   This function requires the libraries ade4 and vegan. If regression-based
#   forward selection is desired, as in the default option, packfor (by Stephane 
#   Dray) is needed. A binary version of packfor for PC can be downloaded from 
#   S. Dray's webpage (http://biomserv.univ-lyon1.fr/~dray/). Mac users will 
#   have to download and compile the source file for their OS level.
#
#   Y:   an object of matrix type containing response data 
#        ready for the analysis (i.e. pretransformed if necessary)
#   space: a matrix containing the geographical (x or x-y) coordinates 
#        of the sites.
#   truncd is a user-provided truncation distance (for cases where the
#        otherwise automatically computed largest value of minimum spanning tree
#        is not appropriate).
#   PCNM: an optional matrix containing precomputed spatial variables like 
#        PCNM, MEM or polynomials.
#   method: specifies the method of selection of PCNM variables: 
#     "none": no selection. All the spatial variables are used. Sometimes
#        useful when the user provides precomputed spatial variables.
#     "fwd": regression-based forward selection (fastest method,
#        default, double stopping criterion [alpha level and adjusted R2], 
#        recommended choice)
#     "AIC": Akaike information criterion (very liberal)
#     "all": test of each PCNM variable in turn, holding all other variables
#        constant (rather conservative)
#     "holm": same as "all" but with Holm correction (very conservative). 
#        Options "all" and "holm" are computationally more intensive
#        and therefore take more time than the other methods of selection.
#     BEWARE: if the rda test involving all the PCNM variables is
#        nonsignificant at the preselected alpha level, the procedure 
#        is stopped, because of the high risk of type I error associated
#        with an a priori selection of explanatory variables. Users of the
#        PCNM method MUST run an overall test prior to ANY type of variable
#        selection. The overall test has a correct type I error
#        [Borcard D. & P. Legendre, Ecological modelling 153(2002):51-68].        
#        Furthermore, a selection of variables should never produce a model
#        with an adjusted R-square higher than the value obtained for the
#        global model containing all the variables. If it happens, a warning
#        is issued. If the difference is important, the user should try a more
#        conservative selection method.
#   alpha: level of significance of the tests (default = 0.05).
#   rangexy: if set to TRUE, rescales the spatial coordinates to a minimum of
#       0 on both axes and a maximum of 10 on the axis wih the largest range,
#       without altering the ratio between the x and y coordinates. If only
#       one coordinate is provided, rescales it into a range of [0;10].
#       Useful when the x-y coordinates are provided in ridiculously large or
#       small values.
#   detrend:  if TRUE (default), prior ot the PCNM analysis, the function tests 
#       if a linear trend is present and significant at the level given by 
#       "alpha", and, if yes, runs a linear detrending of the response data. 
#       It is useless to waste sine-shaped PCNM variables to model linear trends.
#       
#   The function provides spatial plots of the several first canonical 
#        axes (thanks to Guillaume Blanchet for improving the plots!)
#        and a list of results containing the PCNM variables, RDA results 
#        and a global test of the RDA.
#   Several summary diagnostics are displayed on-screen.
#
#   Call (examples): 
#
#       A. Usual case with defaut settings (works well in most situations):
#          name_of_object <- quickPCNM(Y,space)
#       B. Case with AIC selection:
#          name_of_object <- quickPCNM(Y,space,method="AIC")
#       C. Same as A but with a modified alpha level and rescaling of coordinates:
#          name_of_object <- quickPCNM(Y,space,alpha=0.01,rangexy=TRUE)
#       D. Case with precomputed PCNM variables, Holm selection and a modified 
#          alpha level:
#          name_of_object <- quickPCNM(Y,space,method="holm",PCNM=PCNM,alpha=0.01)
#       In all cases it is safer to respect the order of the arguments.
#
#   If you want to run the function several times to compare results, don't
#   forget to ask for a new graphics windows. Otherwise the active window
#   created during the first run will be overwritten.
#   When the run is completed type 'summary(name_of_object)' to get the RDA 
#   results.
#
#                                        Daniel Borcard
#                                        Universite de Montreal
#                                        May 2007 - November 2008

require(ade4)
require(vegan)

Y <- as.matrix(Y)
space <- as.matrix(space)
n <- nrow(Y)

# ----------------------------------------------------------------------------

if(rangexy==TRUE) {

if(ncol(space)==1){
space <- (space-min(space))/((max(space)-min(space))*0.1)
                  }
else{
mini <- apply(space,2,"min")
maxi <- apply(space,2,"max")
xy.trans <- sweep(space,2,mini)
range.max <- 0.1*(max((maxi[1]-mini[1]),(maxi[2]-mini[2])))
space <- as.matrix(xy.trans/range.max)
    }
                }

# ----------------------------------------------------------------------------

if (is.null(PCNM)) {

### 1. Building the PCNM variables

## Computation of a matrix of Euclidean distances among sites
dist.d1 <- dist(space)

## Use of user-defined truncation value for distance matrix
if(truncd > 0) {
  dmin = truncd }

else {
## If no truncation distance is given by the user, search for truncation value 
#  (largest value of minimum spanning tree) using vegan:
spanning <- spantree(dist.d1)
dmin <- max(spanning$dist)
     }

## Truncation of distance matrix - old method with loops
# mat.dist <- as.matrix(dist.d1)
# temp <- matrix(0,nrow(space),nrow(space))
# for(i in 1:nrow(space)){
#     for(j in 1:nrow(space)){
#         temp[i,j] <- ifelse(mat.dist[i,j]>dmin,dmin*4,mat.dist[i,j])
#                            }
#                        }
# dist.trunc <- as.dist(temp)

# Shorter and faster truncation method - not yet fully tested, feel free to try it!
dist.d1[dist.d1 > dmin] <- 4*dmin
dist.trunc <- dist.d1

## Computation of PCNM (with due care to negative eigenvalues)
wa.old <- options(warn = -1)
    on.exit(options(wa.old)) # no warnings for negative eigenvalues
pcoord <- cmdscale(dist.trunc,k=nrow(space)-1,eig=T)
ev <- pcoord$eig

## Counting the positive eigenvalues and building the PCNM matrix
nb.ev <- length(which(ev > 0.0000001))
PCNM <-as.data.frame(pcoord$points[1:nrow(space),1:nb.ev])
PCNM <<- PCNM   # compensation for internal loss of object by R
                     }
else {
PCNM <- as.data.frame(PCNM)
dmin <- "implicit in PCNM file"
meanPCNM <- sum(apply(PCNM,2,mean))
if(abs(meanPCNM) > 0.0000000001 ) {
cat("\n------------------------------------------------------------------")
cat("\nWARNING: the user-provided spatial variables are not centred")
cat("\nto zero mean. Are you sure that they are correct?")
cat("\n------------------------------------------------------------------")
                                  }
sumcorPCNM <- sum(cor(PCNM))
if(abs(sumcorPCNM) > ncol(PCNM)+0.0000000001 ) {
cat("\n------------------------------------------------------------------")
cat("\nWARNING: the user-provided spatial variables are not orthogonal")
cat("\nto one another. Are you sure that they are correct?")
cat("\n------------------------------------------------------------------")
cat("\n")
                                    }
nb.ev <- ncol(PCNM)
ev <- apply(PCNM^2,2,sum)
     }

# ----------------------------------------------------------------------------

### RDA "response matrix x PCNM"

## Preliminary step: linear detrending of response data if trend is significant
if(detrend==TRUE){
   trace.correl <-sum(diag(cor(Y,space)))
   if(abs(trace.correl)<1e-12){
     Y.det <<- Y              
     temp2.test<-matrix(rep(1,5),1)
                              }
   else{  
   temp2 <- rda(Y,space)
   temp2.test <- anova.cca(temp2,alpha=alpha,step=100,perm.pmax=1000)
   if(temp2.test[1,5] <= alpha) {
      Y.det <<- resid(lm(Y~space))  # compensation for internal loss of object by R
       }                 
   else {
      Y.det <<- Y   # compensation for internal loss of object by R
      temp2.test<-matrix(rep(1,5),1)
        }
                 }
                              }
else {
   Y.det <<- Y   # compensation for internal loss of object by R
     }

## RDA with all PCNM variables(complete model)

mod1 <- rda(Y.det~.,data=PCNM)
global.test <- anova.cca(mod1,alpha=alpha,step=100,perm.max=1000)
mod1.sum <- summary(mod1,scaling=1)
R2glob <- mod1.sum$constr.chi/mod1.sum$tot.chi
R2glob.a <- 1-((n-1)/(n-global.test[1,1]-1))*(1-R2glob)
thresh <- R2glob.a+0.001

if(global.test[1,5] >= alpha) {
   cat("\n------------------------------------------------------------------")
   cat ("\n*** Procedure stopped ***")
   cat ("\np-value of global test: ",global.test[1,5])
   cat ("\nNo significant spatial structure detected by global PCNM analysis")
   cat ("\nSelection of PCNM variables would lead to spurious model")
   cat("\n------------------------------------------------------------------","\n")
     }

# ----------------------------------------------------------------------------

else {

METHODS <- c("none","AIC", "fwd", "all","holm")
    method <- match.arg(method, METHODS)

if(method == "none"){ 
mod <- mod1
mod.test <- global.test
mod.sum <- mod1.sum
R2 <- R2glob
R2adj <- R2glob.a
vars.sign <- c(1:ncol(PCNM))
nb.sig.ev <- length(vars.sign)
                    }

else{

if(method == "AIC"){   # 1.1 open  AIC

## Selection of significant PCNM variables using AIC
##     and RDA on selected PCNM variables
mod0 <- rda(Y.det~1,data=PCNM)
mod <- step(mod0,scope=formula(mod1))
mod.sum <- summary(mod,scaling=1)
nb.sig.ev <- nrow(mod.sum$biplot)
cat("\nAIC value for selected model: ",extractAIC(mod)[2],"\n")
                   }   # 1.1 close  AIC
else {                 # 1.2 open    everything but AIC

  if(method == "fwd"){   # 1.2.1 open  fwd

  ## Regression-based forward selection of PCNM variables and RDA on significant
  ##     PCNM variables

  require(packfor)  #  packfor version 0.0-7 or later (to allow adjR2thresh)
   
  # If there is only one response variable that is normally distributed, save
  # time by replacing permutation tests by parametric tests. Otherwise and for
  # a multivariate response matrix, forward selection with the adjusted R2 as 
  # aditional stopping criterion.
    if(ncol(Y)==1){
    norm.test <- shapiro.test(Y.det)
       if(norm.test[2] > 0.05) {
       cat("\n------------------------------------------------------------------")
       cat("\nOnly one, normally distributed response variable found. Parametric")
       cat("\nforward selection on standardized response variable has been run.")
       cat("\n------------------------------------------------------------------","\n")
       Y.det <-scale(Y.det)
       fwd.sel <- forward.sel.par(Y.det,PCNM,alpha=alpha,adjRsqSup=R2glob.a)
                               } else {
       cat("\n------------------------------------------------------------------")
       cat("\nThe only response variable is not normally distributed.")
       cat("\nPermutational forward selection has been run.")
       cat("\n------------------------------------------------------------------","\n")
       fwd.sel <- forward.sel(Y.det,PCNM,alpha=alpha,adjR2thresh = R2glob.a)
                                      }
                  } else {
  fwd.sel <- forward.sel(Y.det,PCNM,alpha=alpha,adjR2thresh = R2glob.a)
                         }
  nb.sig.ev <- nrow(fwd.sel)
  vars.sign <- sort(fwd.sel[,2])}   # 1.2.1 close  fwd

    else {                 # 1.2.2 open   all + holm (common part)
    ## Test of each PCNM variable in turn
    PCNM.pval <- c(rep(0,nb.ev))

      if(method == "all"){               # 1.3 op  all
      for(i in 1:nb.ev){
      PCNM.rda <- rda(Y.det,PCNM[,i],PCNM[,-i],data=PCNM)
      PCNM.test <- anova.cca(PCNM.rda, alpha=alpha, step=100, perm.max=1000,model="reduced")
      PCNM.pval[i] <- PCNM.test[1,5]    }
      vars.sign <- which(PCNM.pval<=alpha)
      nb.sig.ev <- length(vars.sign) }     # 1.3 cl  all
      else{                              # 1.4 op  holm
        # To get enough permutations to allow for Holm correction:
        nperm <- 40*nb.ev
        for(i in 1:nb.ev){
        PCNM.rda <- rda(Y.det,PCNM[,i],PCNM[,-i],data=PCNM)
        PCNM.test <- anova.cca(PCNM.rda, alpha=alpha, step=nperm, perm.max=nperm,model="reduced")
        PCNM.pval[i] <- PCNM.test[1,5]   }
        # Retain significant PCNM variables with Holm correction
        tab.prob <- matrix(cbind(c(1:nb.ev),PCNM.pval),nb.ev)
        tab.prob.tri <-tab.prob[order(tab.prob[,2]),]
        for (i in 1:nb.ev){
        tab.prob.tri[i,2] <- tab.prob.tri[i,2]*(nb.ev-i+1)
                        }
        tab.prob.tri <- tab.prob.tri[order(tab.prob.tri[,2]),]
        # Select variables with p-values <= alpha, with an exception if no such
        # variable has been found.
        if(tab.prob.tri[1,2]>alpha){
    cat("\n------------------------------------------------------------------")
    cat("\nWARNING: despite a significant global test, no single significant")
    cat("\nPCNM variable has been found in Holm-corrected individual tests.")
    cat("\nPlease try with another, less conservative selection procedure.")
    cat("\n------------------------------------------------------------------")
        nb.sig.ev <- 0
        err.no <- options("show.error.messages" = FALSE)
        on.exit(options(err.no))

        stop()

        } else{
        vars.sign <- sort(tab.prob.tri[1:length(which(tab.prob.tri[,2]<=alpha)),1])
        nb.sig.ev <- length(vars.sign)
              }
            }               # 1.4 cl  Holm
         }                  # 1.2.2 cl    all+holm (common part)

#if(method != "AIC" && nb.sig.ev > 0) {

PCNMred <<- as.data.frame(PCNM[,c(vars.sign)])
mod <- rda(Y.det~.,data=PCNMred)
mod.sum <- summary(mod,scaling=1)}

mod.test <- anova.cca(mod, alpha=alpha, step=100,perm.max=1000)

R2 <- mod.sum$constr.chi/mod.sum$tot.chi
R2adj <- 1-((n-1)/(n-mod.test[1,1]-1))*(1-R2)

   }                        # close all choices (incl. no selection)

# Warning if the adjusted R-square of minimal model is greater than the
# adjusted R-square of the global model with all PCNM variables.

if(R2adj > (R2glob.a+0.05*R2glob.a)){
cat("\n------------------------------------------------------------------")
cat("\nWARNING: the adjusted R-square of the reduced model,",round(R2adj,4))
cat("\nexceeds the adjusted R-square of the global model,",round(R2glob.a,4))
cat("\nby more that 5%.This means that the selection has been overly liberal.")
cat("\nChoose another, mode conservative method (see explanations).")
cat("\n------------------------------------------------------------------")
                                   }
else if(R2adj > (R2glob.a) & R2adj <= (R2glob.a+0.05*R2glob.a)){
cat("\n------------------------------------------------------------------")
cat("\nWARNING: the adjusted R-square of the reduced model,",round(R2adj,4))
cat("\nexceeds the adjusted R-square of the global model,",round(R2glob.a,4))
cat("\nby 5% or less. This means that the selection has been a little bit")
cat("\ntoo liberal. This small amount should not harm, however.")
cat("\n------------------------------------------------------------------")
                                   }

# At this point, anova.cca for all axes doesn't seem to work properly: it looses
# the track of one or the other object defined higher (generally PCNMred).
# The workaround is to export PCNMred from the function during the run
# using the trick: PCNMred <<- ...  (Line 343). PCNMred is thus present in the
# main R workspace and should be deleted prior to another quickPCNM run.

 mod.axes.test <- anova.cca(mod,by="axis",step=100,perm.max=1000)
 
## Count the significant axes 
 nb.ax <- 0
 i <- 1
    while(mod.axes.test[i,5]<=alpha){
       nb.ax=nb.ax+1
       i = i+1
                 }

# ----------------------------------------------------------------------------

## Number of axes to draw (arbitrary rule, an alterative to the tests above)
#  if(mod.test[1,1]<=2){
#     nb.ax=mod.test[1,1]} else {
#     nb.ax=ceiling(sqrt(mod.test[1,1]))}

## Plot of significant axes
fitted.scores=scores(mod,display="lc",choices=1:nb.ax)
par(mfrow=c(round(sqrt(nb.ax)),ceiling(sqrt(nb.ax))))

if(ncol(space)==2){
   for(i in 1:nb.ax){
   s.value(space,fitted.scores[,i],addaxes=FALSE,include.origin=FALSE,sub=paste("Axis ",i),csub=1.5)
                    }
} else {
   for(i in 1:nb.ax){
   plot(space,fitted.scores[,i],type="l",ylab="Fitted site scores")
                    }
       }

# ----------------------------------------------------------------------------

## Screen output
if(detrend==TRUE){
   if(temp2.test[1,5] <= alpha) {
      cat("\n-------------------------------------------------------")
      cat("\nA significant linear trend has been found in the response data.")
      cat("\nThe data have been detrended prior to PCNM analysis.")
                                }
    else{
      cat("\n-------------------------------------------------------")
      cat("\nNo significant linear trend has been found in the response data.")
      cat("\nThe data have NOT been detrended prior to PCNM analysis.")     
        }
                 }


cat("\n-------------------------------------------------------")
cat("\nThe truncation value used for PCNM building is",dmin)
cat("\nThere are ",nb.ev," positive PCNM eigenvalues","\n")
cat("Adjusted R2 of global model = ",round(R2glob.a,4),"\n")
if(method != "none") {
if(nb.sig.ev==1){
cat(nb.sig.ev," PCNM variable has been selected","\n")}
else{
cat(nb.sig.ev," PCNM variables have been selected","\n")}
cat("R2 of minimum model = ",round(R2,4),"                    ","\n")
cat("Adjusted R2 of minimum model = ",round(R2adj,4),"        ","\n")
                     }
# cat("The minimum model has ",nb.ax,"significant canonical axes","\n")
cat("---------------------------------------------------------")
cat("\n")
## Extraction of results to be returned
#table <- list(PCNM,mod.sum,mod.test,mod.axes.test)
#names(table) <- c("PCNM","RDA","RDA_test","Ax_tests")
if(method == "none") {
   table <- list(mod.sum,mod.test)
   names(table) <- c("RDA","RDA_test")
                     }
else{

   if(method == "AIC"){
   table <- list(PCNM,ev[1:nb.ev],mod.sum,mod.test)
   names(table) <- c("PCNM","PCNM_eigenvalues","RDA","RDA_test")
                      }
   else{
   if(method == "fwd"){
   table <- list(PCNM,ev[1:nb.ev],fwd.sel,PCNMred,mod.sum,mod.test)
   names(table) <- c("PCNM","PCNM_eigenvalues","fwd.sel","PCNM_reduced_model","RDA","RDA_test")
                      }
   else{
   table <- list(PCNM,ev[1:nb.ev],PCNMred,mod.sum,mod.test)
   names(table) <- c("PCNM","PCNM_eigenvalues","PCNM_reduced_model","RDA","RDA_test")
                      }
                      }
    }
return(table)
}
}

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

forward.sel.par <- function(Y, X, alpha = 0.05, K = nrow(X)-1, R2thresh = 0.99, R2more = 0.001, adjRsqSup = 0.99, Yscale = FALSE, verbose=TRUE)
#
# Parametric forward selection of explanatory variables in regression and RDA.
# Y is the response, X is the table of explanatory variables.
#
# If Y is univariate, this function implements FS in regression.
# If Y is multivariate, this function implements FS using the F-test described 
# by Miller and Farr (1971). This test requires that
#   -- the Y variables be standardized,
#   -- the error in the response variables be normally distributed (to be verified by the user).
#
# This function uses 'simpleRDA2' and 'RsquareAdj' developed for 'varpart' in 'vegan'.
#
#                Pierre Legendre & Guillaume Blanchet, May 2007
#
# Arguments --
#
# Y         Response data matrix with n rows and m columns containing quantitative variables.
# X         Explanatory data matrix with n rows and p columns containing quantitative variables.
# alpha     Significance level. Stop the forward selection procedure if the p-value of a variable is higher than alpha. The default is 0.05.
# K         Maximum number of variables to be selected. The default is one minus the number of rows.
# R2thresh  Stop the forward selection procedure if the R-square of the model exceeds the stated value. This parameter can vary from 0.001 to 1.
# R2more    Stop the forward selection procedure if the difference in model R-square with the previous step is lower than R2more. The default setting is 0.001.
# adjRsqSup Stop the forward selection procedure if the adjusted R-square of the model exceeds the stated value. This parameter can take any value (positive or negative) smaller than 1.
# Yscale    Standardize the variables in table Y to variance 1. The default setting is FALSE. The setting is automatically changed to TRUE if Y contains more than one variable. This is a validity condition for the parametric test of significance (Miller and Farr 1971).
#
# Reference:
# Miller, J. K., and S. D. Farr. 1971. Bimultivariate redundancy: a comprehensive measure of 
#    interbattery relationship. Multivariate Behavioral Research 6: 313-324.

{
library(vegan)
Y = as.matrix(Y)
X = apply(as.matrix(X),2,scale,center=TRUE,scale=TRUE)
var.names = colnames(as.data.frame(X))
n = nrow(X)
m = ncol(X)
if(nrow(Y) != n) stop("Numbers of rows not the same in Y and X")
p = ncol(Y)
if(p > 1) {
   Yscale = TRUE
   if(verbose) cat("The variables in response matrix Y have been standardized",'\n')
   }
Y = apply(Y,2,scale,center=TRUE,scale=Yscale)
SS.Y = sum(Y^2)

X.out = c(1:m)

# Find the first variable X to include in the model
R2prev = 0
R2cum = 0
for(j in 1:m) {
   toto = simpleRDA2(Y,X[,j],SS.Y)
   if(toto$Rsquare > R2cum) {
      R2cum = toto$Rsquare
      no.sup = j
      }
   }
mm = 1
FP = FPval(R2cum,R2prev,n,mm,p)
if(FP$pval <= alpha) {
   adjRsq = RsquareAdj(R2cum,n,mm)
   res1 = var.names[no.sup]
   res2 = no.sup
   res3 = R2cum
   res4 = R2cum
   res5 = adjRsq
   res6 = FP$Fstat
   res7 = FP$pval
   X.out[no.sup] = 0
   delta = R2cum
   } else { stop("Procedure stopped (alpha criterion): pvalue for variable ",no.sup," is ",FP$pval) }

# Add variables X to the model
while((FP$pval <= alpha) & (mm <= K) & (R2cum <= R2thresh) & (delta >= R2more) & (adjRsq <= adjRsqSup)) {
   mm = mm+1
   R2prev = R2cum
   R2cum = 0
   for(j in 1:m) {
      if(X.out[j] != 0) {
         toto = simpleRDA2(Y,X[,c(res2,j)],SS.Y)
         if(toto$Rsquare > R2cum) {
            R2cum = toto$Rsquare
            no.sup = j
            }
         }
      }
   FP = FPval(R2cum,R2prev,n,mm,p)
   delta = R2cum-R2prev
   adjRsq = RsquareAdj(R2cum,n,mm)
   res1 = c(res1,var.names[no.sup])
   res2 = c(res2,no.sup)
   res3 = c(res3,delta)
   res4 = c(res4,R2cum)
   res5 = c(res5,adjRsq)
   res6 = c(res6,FP$Fstat)
   res7 = c(res7,FP$pval)
   X.out[no.sup] = 0
   }
if(verbose) {
	if(FP$pval > alpha)  cat("Procedure stopped (alpha criterion): pvalue for variable ",no.sup," is ",FP$pval,'\n') 
	if(mm > K)           cat("Procedure stopped (K criterion): mm = ",mm," is larger than ",K," after including variable ",no.sup,'\n') 
	if(R2cum > R2thresh) cat("Procedure stopped (R2thresh criterion): R2cum for variable ",no.sup," is ",R2cum,'\n') 
	if(delta < R2more)   cat("Procedure stopped (R2more criterion): delta for variable ",no.sup," is ",delta,'\n') 
	if(adjRsq>adjRsqSup) cat("Procedure stopped (adjRsqSup criterion): adjRsq for variable ",no.sup," is ",adjRsq,'\n') 
	}

res=data.frame(res1,res2,res3,res4,res5,res6,res7)
colnames(res) = c("variable","order","R2","R2cum","AdjRsq","F","pval")
if((FP$pval > alpha) | (mm > K) | (R2cum > R2thresh) | (delta < R2more) | (adjRsq > adjRsqSup))  res = res[1:(mm-1),]

return(res)
}

FPval <- function(R2cum,R2prev,n,mm,p)
# Compute the partial F and p-value after adding a single explanatory variable to the model.
# In FS, the number of df of the numerator of F is always 1. See Sokal & Rohlf 1995, eq 16.14.
# 
# The amendment, based on Miller and Farr (1971), consists in multiplying the numerator and  
# denominator df by 'p', the number of variables in Y, when computing the p-value.
#
#                Pierre Legendre, May 2007
{
df2 = (n-1-mm)
Fstat = ((R2cum-R2prev)*df2) / (1-R2cum)
pval = pf(Fstat,1*p,df2*p,lower.tail=FALSE)
return(list(Fstat=Fstat,pval=pval))
}