manova.2way.unbalanced <-
	function(Y, A, B, nperm=999, model="direct")
#
# 2-way crossed-factor manova for univariate or multivariate response 
# data, balanced or unbalanced designs, with permutation tests. 
# This is a Model I manova, meaning that the two factors are fixed.
#
# For unbalanced designs, this function computes type III sums-of-squares. 
# They are equal to type I sums of squares in the case of balanced designs. 
# The function can take either a single variable or a whole data table as 
# the response Y. The computation method by RDA was described by Legendre 
# and Anderson (1999) for balanced designs.
#
# This function requires replicates to test the interaction. 
# It cannot compute a two-way manova without replication.
#
# Arguments -
#
#   Y = vector or matrix of response variable(s)
#   A = factor A
#   B = factor B
#   nperm = number of permutations. Default value: nperm=999
#   model = c("reduced", "direct", "full")
#         Permutation model in vegan's 'anova.cca' function.
#         Default value: model="direct". Following Anderson and Legendre (1999),
#         permutation of the raw data is adequate for anova since there are no 
#         outlier values in the factors.
#
# References -
# 
# Anderson, M. J. and P. Legendre. 1999. An empirical comparison of permutation
# methods for tests of partial regression coefficients in a linear model. 
# Journal of Statistical Computation and Simulation 62: 271-303.
#
# Legendre, P. & M. J. Anderson. 1999. Distance-based redundancy analysis: 
# testing multispecies responses in multifactorial ecological experiments. 
# Ecological Monographs 69 (1): 1-24.
#
#                                                  Pierre Legendre, October 2008
#
################################################################################
{
library(vegan)
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
out.rda <- rda(Y, model.A, cbind(model.B, model.int))
int.out <- anova(out.rda, perm.max=(nperm+1), step=(nperm+1), model=model)
out.A <- int.out[1,]

# Test of factor B
out.rda <- rda(Y, model.B, cbind(model.A, model.int))
int.out <- anova(out.rda, perm.max=(nperm+1), step=(nperm+1), model=model)
out.B <- int.out[1,]

# Test of interaction
out.rda <- rda(Y, model.int, cbind(model.A, model.B))
int.out <- anova(out.rda, perm.max=(nperm+1), step=(nperm+1), model=model)
out.AB <- int.out[1,]
out.res <- int.out[2,]

# Parametric probabilities (for univariate response data only)
if(p == 1) {
	P.param <- pf(out.A[1,3],out.A[1,1],out.res[1,1],lower.tail=FALSE)
	P.param <- c(P.param, pf(out.B[1,3],out.B[1,1],out.res[1,1],lower.tail=FALSE))
	P.param <- c(P.param, pf(out.AB[1,3],out.AB[1,1],out.res[1,1],lower.tail=FALSE))
	P.param <- c(P.param, NA)
	}

# Table of results
out <- rbind(out.A, out.B, out.AB, out.res)
out$Var <- out$Var*(n-1)/out$Df
rownames(out) <- c("Factor A","Factor B","Interaction","Residuals")
if(p == 1) {
	out <- as.data.frame(cbind(out[,1:4], P.param, out[,5]))
	colnames(out) <- c("Df","MeanSq","F-stat","N.perm","P.param","Pr(>F)")
	} else {
	colnames(out) <- c("Df","MeanSq","F-stat","N.perm","Pr(>F)")
	}

out <- list(Anova_unbalanced=out)
class(out) <- "anova.unbalanced"
return(out)
}
