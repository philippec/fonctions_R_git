t.perm <- function(vec1, vec2, nperm=999, tail=2, silent=FALSE)
#
# This function computes a permutation test of comparison of the means
# of two vectors corresponding to independent samples.
#
# The formula for the weighted mean of the variances is used in the
# calculation of t. So, the t-statistic here may differ slightly from that
# of a t-test with Welch correction, where that correction is not used. 
# The Welch correction for t-tests is only available for parametric tests.
# There is no known permutational equivalent.
#
# The permutational t-test requires equality or near-equality of the 
# within-group variances. It does not require normality of the distributions.
#
# Parametres of the function
#
#    vec1, vec2: the two vectors to be compared
#    nperm = number of permutations (default value: 999)
#    tail = -1 : one-tailed test in the lower tail
#    tail =  1 : one-tailed test in the upper tail
#    tail =  2 : two-tailed test (default value)
#    silent = FALSE: calculation results are printed to the R console.
#           = TRUE : calculation results are not printed to the R console
#                    (useful for simulations).
#
# Values returned
# 
#    t.ref : reference value of the t-statistic
#    p.param : parametric p-value
#    p.perm : permutational p-value
#    nperm : number of permutations
#    perm.t : list of the t statistics, starting with the reference value,
#             followed by all values obtained under permutations.
#
# Example: t-test between two vectors of random normal deviates
#
#     group1 <- rnorm(50,0,10)              # Group 1 has a mean of 0
#     group2 <- rnorm(30,5,10)              # Group 2 has a mean of 5
#     res = t.perm(group1, group2, 999, -1) # One tailed test in the lower tail
#
#   Compare with results of the parametric t-test:
# 
#     t.test(group1,group2,var.equal=TRUE,alternative="less") 
#
# Pedagogic objectives --
# 
# - Show the  structure of a function in R
# - Show hot to program a 'for' loop: "for(i in 1:nperm)"
# - Show how to print out results to the R console: function "cat"
# - Show how to program a permutation test, one-tailed or two-tailed
#
#                          Pierre Legendre, November 2005
{
# Find the lengths of the two vectors
n1 <- length(vec1)
n2 <- length(vec2)
n  <- n1+n2

# Combine the two vectors into a single vector 'vec'
vec <- c(vec1,vec2)

# Compute the t-statistic for the unpermuted data
t.ref <- t.stat(n1,n2,vec1,vec2)
perm.t <- t.ref$stat   # Start the vector containing the list of t statistics

# Compute the parametric p-value
if(tail == -1) p.param <- pt(t.ref$stat,(n-2),lower.tail=TRUE)
if(tail ==  1) p.param <- pt(t.ref$stat,(n-2),lower.tail=FALSE)
if(tail ==  2) p.param <- pt(abs(t.ref$stat),(n-2),lower.tail=FALSE)*2

# Print these first results
if(!silent) cat('\nt-test comparing two group means','\n')
if(!silent) cat('Group sizes:',n1,n2,'\n')
if(!silent) cat('Group means:',t.ref$moy1,t.ref$moy2,'\n')
if(!silent) cat('Group variances:',t.ref$var1,t.ref$var2,'\n')
if(!silent) cat('t =',t.ref$stat,'\n')
if(!silent) cat('d.f. =',(n-2),'\n')
if(!silent) cat('Prob (parametric) =',p.param,'\n')

# Perform the permutation test
nPGE <- 1
for(i in 1:nperm)
   {
   vec.perm  <-  sample(vec,n)
   vec1.perm <- vec.perm[1:n1]
   vec2.perm <- vec.perm[(n1+1):n]
   t.perm <- t.stat(n1,n2,vec1.perm,vec2.perm)
   perm.t <- c(perm.t, t.perm$stat)
   
   if(tail == -1) if(t.perm$stat <= t.ref$stat) nPGE <- nPGE+1
   if(tail ==  1) if(t.perm$stat >= t.ref$stat) nPGE <- nPGE+1
   if(tail ==  2) if( abs(t.perm$stat) >= abs(t.ref$stat) ) nPGE <- nPGE+1
   }

# Compute and print the permutational p-value
P <- nPGE/(nperm+1)
if(!silent) cat('Prob (',nperm,'permutations) =',formatC(P,digits=5,width=7,format="f"),'\n')
#
return(list(t.ref=t.ref$stat, p.param=p.param, p.perm=P, nperm=nperm, perm.t=perm.t))
}

t.stat <- function(n1,n2,vec1,vec2)
# Compute the t-statistic 
{
moy1 <- mean(vec1)
moy2 <- mean(vec2)
var1 <- var(vec1)
var2 <- var(vec2)
var.wm <- ( (n1-1)*var1 + (n2-1)*var2 ) / (n1+n2-2)
t <- ( moy1-moy2 ) / sqrt(var.wm * ((1/n1) + (1/n2)) )
return(list(moy1=moy1,moy2=moy2,var1=var1,var2=var2,stat=t))
}
