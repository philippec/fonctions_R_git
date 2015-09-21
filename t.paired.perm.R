t.paired.perm <- function(vec1, vec2, nperm=999, alternative="two.sided", silent=FALSE)
#
# This function computes a permutation test of comparison of the means
# of two paired vectors (related samples).
# For each object, permutations are restricted to the two related observations.
#
# Parametres of the function
#
#    vec1, vec2: the two vectors to be compared
#    nperm = number of permutations (default value: 999)
#    alternative = c("two.sided", "less", "greater"). Dafault value: "two.sided"
#    silent = FALSE: calculation results are printed to the R console.
#           = TRUE : calculation results are not printed to the R console (for simulations).
#
# Values returned
# 
#    t.ref : reference value of the t-statistic
#    p.param : parametric p-value
#    p.perm : permutational p-value
#    nperm : number of permutations
#
# Example: Deer leg length data from Zar (1999, p. 162)
#
# deer = matrix(c(142,140,144,144,142,146,149,150,142,148,138,136,147,139,143,141,143,145,136,146), 10, 2)
# rownames(deer) = c('Deer.1','Deer.2','Deer.3','Deer.4','Deer.5','Deer.6','Deer.7','Deer.8','Deer.9','Deer.10')
# colnames(deer) = c('Hind.leg', 'Fore.leg')
#
# res = t.paired.perm(deer[,1], deer[,2])   # Two-tailed test by default
#
# Compare results to:  res2 = t.test(deer[,1], deer[,2], paired=TRUE) 
#
#                          Pierre Legendre, November 2009
#                          Guillaume Blanchet, September 2015 (faster permutation code)
{
n1 <- length(vec1)
n2 <- length(vec2)
if(n1 != n2) stop("The two vectors have different lengths. They cannot be paired.")

tail <- match.arg(alternative, c("two.sided", "less", "greater"))

res = t.test(vec1, vec2, paired=TRUE, alternative=tail)
t.ref =  res$statistic

# Print these first results
if(!silent) cat('\nt-test comparing the means of two related samples','\n','\n')
if(!silent) cat('Number of objects:',n1,'\n')
if(!silent) cat('Mean of the differences:',res$estimate,'\n')
if(!silent) cat('t statistic (paired observations):',t.ref,'\n')
if(!silent) cat('95 percent confidence interval of t:',res$conf.int,'\n')
if(!silent) cat('Degrees of freedom:',res$parameter,'\n')
if(!silent) cat('Alternative hypothesis:',tail,'\n')
if(!silent) cat('Prob (parametric):',res$p.value,'\n')

# Perform the permutation test
# Permutations are restricted to the two related observations for each object.
nPGE <- 1

for(i in 1:nperm)
	{
## New permutation code, GB
	mat <- cbind(vec1,vec2)
	topermute <- rbinom(n1,1,0.5)
	mat[topermute==1,] <- mat[topermute==1,2:1]
## End new code, GB
	
	res.perm = t.test(mat[,1], mat[,2], paired=TRUE, alternative=tail)
	t.perm = res.perm$statistic
	
	if(tail == "two.sided") if( abs(t.perm) >= abs(t.ref) ) nPGE <- nPGE+1
	if(tail == "less")      if(t.perm <= t.ref) nPGE <- nPGE+1
	if(tail == "greater")   if(t.perm >= t.ref) nPGE <- nPGE+1
	}

# Compute and print the permutational p-value
P <- nPGE/(nperm+1)
if(!silent) cat('Prob (',nperm,'permutations):',formatC(P,digits=5,width=7,format="f"),'\n')
#
return(list(t.ref=t.ref, p.param=res$p.value, p.perm=P, nperm=nperm))
}
