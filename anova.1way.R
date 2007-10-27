anova.1way <- function(formula, data = NULL, nperm=0)
#
################################################################################
#
# One-way anova with permutation test.
# The factor may be fixed or random, the design balanced or unbalanced.
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
# Notes   : 1. If 'data' is missing, the model must be specified as  Y~X  
#              where Y is the response and X is the anova factor.
#           2. Make sure X is defined 'as.factors', unless it is binary.
#           3. Permutation of the raw data is adequate for anova.
#
#                                                   Pierre Legendre, August 2007
#
################################################################################
#
## Example modified from Venables and Ripley (2002) p. 165.
#
# F1 = as.factor(c(0,1,0,1,1,1,0,0,0,1,1,0,1,1,0,0,1,0,1,0,1,1,0,0))
# yield <- c(49.5,62.8,46.8,57.0,59.8,58.5,55.5,56.0,62.8,55.8,69.5,55.0,
#            62.0,48.8,45.5,44.2,52.0,51.5,49.8,48.8,57.2,59.0,53.2,56.0)
# mat = as.data.frame(cbind(yield,F1))
# 
## Compute one-way anova with permutation test.
#
# out.1way = anova.1way(yield~F1, data=mat, nperm=999)
#
################################################################################
{
# The reference one-way anova
toto = lm(formula, data)
anova.res = anova(toto)

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
    n = nrow(mf)

k = nrow(anova.res) - 1
note = "One-way anova"

if(nperm > 0) {
   GE = rep(1,k)
   Pperm = c(rep(0,k), NA)
   for(j in 1:nperm)
      {
      Yperm = sample(mf[,1],n)
      toto = lm(Yperm ~ F1)
      anova.per = anova(toto)
      for(i in 1:k) { if(anova.per[i,4] >= anova.res[i,4]) GE[i] = GE[i] + 1 }
      }
   for(i in 1:k) { Pperm[i]=GE[i]/(nperm+1) }
   anova.res  = data.frame(anova.res, Pperm)
   colnames(anova.res) = c("Df", "Sum Sq", "Mean Sq", "F value", "Prob(param)", "Prob(perm)")
   note = "One-way anova with permutation test"
   }

if(nperm <= 0) { 
   return(list(anova.type=note,anova.table=anova.res)) 
   } else {
   return(list(anova.type=note, nperm=nperm, response.var=var.names[1], anova.table=anova.res))
   }
}
