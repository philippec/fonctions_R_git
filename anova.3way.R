anova.3way <- function(formula, data = NULL, nperm=0)
#
################################################################################
#
# Three-way crossed-factor anova (fixed factors, balanced designs) 
# with permutation tests.
#
# Arguments
#
#  formula: a formula specifying the model, as in 'lm' and 'aov'.
#
#  data   : A data frame in which the variables specified in the formula are to  
#           be found. 
# 
#  nperm  : number of permutations. A commonly-used value is 999.
#
# Notes   : 1. If 'data' is missing, the model must be specified as  Y~X1*X2*X3  
#              where Y is the response, X1, X2 and X3 are the three factors.
#           2. Make sure X1, X2 and X3 are defined 'as.factors' unless they are
#              binary.
#           3. Following Anderson and Legendre (1999), permutation of the raw
#              data is adequate for anova since there are no outlier values
#              in the factors.
#
# Reference -
# Anderson, M. J. and P. Legendre. 1999. An empirical comparison of permutation
# methods for tests of partial regression coefficients in a linear model. 
# Journal of Statistical Computation and Simulation 62: 271-303.
#
#                                                   Pierre Legendre, August 2007
#
################################################################################
#
## Example modified from Venables and Ripley (2002) p. 165.
#
# F1 = as.factor(c(0,1,0,1,1,1,0,0,0,1,1,0,1,1,0,0,1,0,1,0,1,1,0,0))
# F2 = as.factor(c(1,1,0,0,0,1,0,1,1,1,0,0,0,1,0,1,1,0,0,1,0,1,1,0))
# F3 = as.factor(c(1,0,0,1,0,1,1,0,0,1,0,1,0,1,1,0,0,0,1,1,1,0,1,0))
# yield <- c(49.5,62.8,46.8,57.0,59.8,58.5,55.5,56.0,62.8,55.8,69.5,55.0,
#            62.0,48.8,45.5,44.2,52.0,51.5,49.8,48.8,57.2,59.0,53.2,56.0)
# mat = as.data.frame(cbind(yield,F1,F2,F3))
# 
## Compute 3-way anova with permutation tests. F1, F2, and F3 are fixed factors.
#
# out.3way = anova.3way(yield~F1*F2*F3, data=mat, nperm=999)
#
################################################################################
{
# The reference 3-way anova (3 fixed factors)
toto = lm(formula, data)
anova.res = anova(toto)
if(anova.res[8,1] == 0) stop("There is no 3-way interaction in this problem.")

# From the formula, find the variables as well as the number of observations 'n'
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action", 
        "offset"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    var.names = colnames(mf)
    Y = mf[,1]
    F1 = as.factor(mf[,2])
    F2 = as.factor(mf[,3])
    F3 = as.factor(mf[,4])
    n = nrow(mf)
    k = nrow(anova.res) - 1

# Check model balance
balance = var(as.vector(table(mf[,2:4])))
if(balance > 0) 
stop("Design unbalanced. This function can only handle balanced designs.")

note = "Three-way anova (three fixed factors)"

if(nperm > 0) {
   GE = rep(1,k)
   Pperm = c(rep(0,k), NA)
   for(j in 1:nperm)
      {
      Yperm = sample(mf[,1],n)
      toto = lm(Yperm ~ F1*F2*F3)
      anova.per = anova(toto)
      for(i in 1:k) { if(anova.per[i,4] >= anova.res[i,4]) GE[i] = GE[i] + 1 }
      }
   for(i in 1:k) { Pperm[i]=GE[i]/(nperm+1) }
   anova.res  = data.frame(anova.res, Pperm)
   colnames(anova.res) = c("Df", "Sum Sq", "Mean Sq", "F value", "Prob(param)", "Prob(perm)")
   note = "Three-way anova (three fixed factors) with permutation tests"
   }

if(nperm <= 0) { 
   return(list(anova.type=note,anova.table=anova.res)) 
   } else {
   return(list(anova.type=note, nperm=nperm, response.var=var.names[1], anova.table=anova.res))
   }
}
