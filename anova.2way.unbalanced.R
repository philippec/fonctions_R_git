anova.2way.unbalanced <-
	function(Y, A, B, nperm=999, model="direct", strata=FALSE, silent=FALSE)
#
# Model I, 2-way crossed-factor multivariate anova (by RDA) for balanced or   
# unbalanced designs, with permutation tests. 
# Model I anova means that the two factors are fixed.
#
# For unbalanced designs, this function computes type III sums-of-squares. 
# They are equal to type I sums of squares in the case of balanced designs. 
# The function can take either a single variable or a whole data table as 
# the response Y. The computation method by RDA was described by Legendre 
# and Anderson (1999) for balanced designs.
#
# Usage -
#
# anova.2way.unbalanced(Y, A, B, nperm=999, model="direct", strata=FALSE, silent=FALSE)
#
# Arguments -
#
#   Y = vector or matrix of response variable(s)
#   A = factor A
#   B = factor B
#   nperm = number of permutations. Default value: nperm=999
#   model = c("reduced", "direct", "full")
#         Permutation model available in vegan's 'anova.cca' function.
#         Default value: model="direct". Following Anderson and Legendre (1999),
#         permutation of the raw data is adequate for anova since there are no 
#         outlier values in the factors.
#   strata = FALSE : permutations are NOT performed within the levels of the other factor.
#   strata = TRUE  : permutations are performed within the levels of the other factor.
#   silent = FALSE : the output message at the beginining of the funciton is printed.
#   silent = TRUE  : the output message is not printed (for example, when the  
#                    function is used in a simulation study).
#
# Details -
#
# Disequilibrium of the design reduces the power of the test of significance.
#
# Use of option "strata": permuting, or not, within the levels of the other factor did 
# not seem to make much difference in the results of the numerical simulations that we 
# conducted.
#
# Value -
#
# An anova table showing the results of the permutation tests of the main factors and the 
# interaction. The parametric p-values are also shown when the analysis implies a single 
# response variable.
#
# References -
# 
# Anderson, M. J. and P. Legendre. 1999. An empirical comparison of permutation
# methods for tests of partial regression coefficients in a linear model. 
# Journal of Statistical Computation and Simulation 62: 271-303.
#
# Legendre, P. & M. J. Anderson. 1999. Distance-based redundancy analysis: 
# testing multispecies responses in multifactorial ecological experiments. 
# Ecological Monographs 69: 1-24.
#
### Examples - 
# 
# A = gl(4, 5)
# B = factor(rep(c(1,1,1,2,2),4))
# 
### Multivariate unbalanced data 
### Random response data -- No significant effect is expected
# Y = matrix(rnorm(40),20,2)
# (res = anova.2way.unbalanced(Y,A,B))
# 
### Univariate unbalanced data
### Random response data -- No significant effect is expected
# y = rnorm(20)
# (res.out = anova.2way.unbalanced(y,A,B))
### The results are identical to those of univariate 2-way unbalanced anova
### See http://mcfromnz.wordpress.com/2011/03/02/anova-type-iiiiii-ss-explained/
### Set the contrasts for each factor in lm()
### Use Anova() from package {car}
# require(car)
# res.lm = lm(y ~ A*B, contrasts=list(A=contr.sum, B=contr.sum))
# (res.out.III = Anova(res.lm, type="III"))
# 
# Author:: Pierre Legendre, Université de Montréal, 2008, 2015
# License: GPL-2
#
################################################################################
{
require(vegan)
if(!silent) cat("Output of this function is adjusted to the changes in vegan 2.2-1\n")
Y = as.matrix(Y)
n <- nrow(Y)
p <- ncol(Y)
A <- as.factor(A)
ncol.A <- nlevels(A)-1
B <- as.factor(B)
ncol.B <- nlevels(B)-1
ncol.int <- ncol.A * ncol.B
model.mat <- model.matrix(~ A*B, contrasts=list(A="contr.helmert", B="contr.helmert"))[,-1]

model.A <- as.matrix(model.mat[,1:ncol.A])
model.B <- as.matrix(model.mat[,(ncol.A+1):(ncol.A+ncol.B)])
model.int <- as.matrix(model.mat[,(ncol.A+ncol.B+1):ncol(model.mat)])

# Test of factor A
out.rda1 <- rda(Y, model.A, cbind(model.B, model.int))
# if(strata) { bl<-B } else { bl=NULL }
if(strata) { 
	int.out1 <- anova(out.rda1, permutations=how(blocks=B, nperm=nperm), model=model)
	} else {
	int.out1 <- anova(out.rda1, permutations=how(nperm=nperm), model=model)
	}
out.A <- int.out1[1,]

# Test of factor B
out.rda2 <- rda(Y, model.B, cbind(model.A, model.int))
# if(strata) { bl<-A } else { bl=NULL }
if(strata) { 
	int.out2 <- anova(out.rda2, permutations=how(blocks=A, nperm=nperm), model=model)
	} else {
	int.out2 <- anova(out.rda2, permutations=how(nperm=nperm), model=model)
	}
out.B <- int.out2[1,]

# Test of interaction
out.rda3 <- rda(Y, model.int, cbind(model.A, model.B))
int.out3 <- anova(out.rda3, permutations=how(nperm=nperm), model=model)
out.AB <- int.out3[1,]
out.res <- int.out3[2,]

# Parametric probabilities (for univariate response data only)
if(p == 1) {
	P.param <- pf(out.A[1,3],out.A[1,1],out.res[1,1],lower.tail=FALSE)
	P.param <- c(P.param, pf(out.B[1,3],out.B[1,1],out.res[1,1],lower.tail=FALSE))
	P.param <- c(P.param, pf(out.AB[1,3],out.AB[1,1],out.res[1,1],lower.tail=FALSE))
	P.param <- c(P.param, NA)
	}

# Table of results
out <- rbind(out.A, out.B, out.AB, out.res)
# SumSq <- out$Variance*(n-1)           # Sums of squares, as in Anova output of {car}
# out$Var <- out$Variance*(n-1)/out$Df  # Mean squares
rownames(out) <- c("Factor A","Factor B","Interaction","Residuals")
if(p == 1) {
#	out <- as.data.frame(cbind(out[,1:3], P.param, out[,4]))
	out <- as.data.frame(cbind(out[,1], round(out[,2:3],4), round(P.param,4), out[,4]))
	colnames(out) <- c("Df","MeanSq","F-stat","P.param","Pr(>F)")
#	} else {
#	colnames(out) <- c("Df","MeanSq","F-stat","Pr(>F)")
	}

out <- list(Anova_unbalanced=out)
class(out) <- "anova.unbalanced"
return(out)
}