corPerm1 <- function(x,y,nperm)
#
# En - This function computes a two-tailed permutation test of a Pearson 
# correlation coefficient between two data vectors. The test statistic is r.
#
# Fr - Cette fonction realise, par permutation, un test bilateral du coefficient 
# de correlation de Pearson entre deux vecteurs. La statistique du test est r.
#
# Parameters of the function:
#    x, y: the two data vectors
#    nperm = number of permutations
#
# Example: test the correlation between two vectors of random numbers N(0,1)
#   x <- rnorm(50,0,1)
#   y <- rnorm(50,0,1)
#   cor.out = corPerm1(x, y, 999)
#   cor.out
# Compare the result to:  cor.test(x,y)
#
# Pedagogic objectives --
# - illustrate the structure of a function in the R language
# - show how to program a "for" loop (command: "for(i in 1:nperm)")
# - show how to write out the results using "cat"
# - show how to program a permutation test
#
#                          Pierre Legendre, October 2005
{
x <- as.matrix(x)
y <- as.matrix(y)
n <- nrow(x)
r.ref <- cor(x,y)

nGT <- 1
for(i in 1:nperm)
   {
   y.perm  <-  sample(y,n)
   r.perm <- cor(x,y.perm)
   if( abs(r.perm) >= abs(r.ref) ) nGT <- nGT+1
   }
P <- nGT/(nperm+1)
cat('\nPearson correlation (two-tailed test)','\n')
cat('r =',r.ref,'\n')
cat('Prob(',nperm,'permutations) =',P,'\n','\n')
return(list( Correlation=r.ref, No.perm=nperm, P.perm=P ))
}


corPerm2 <- function(x,y,nperm)
#
# En - This function computes a two-tailed permutation test of a Pearson 
# correlation coefficient between two data vectors. The test statistic is t.
#
# Fr - Cette fonction realise, par permutation, un test bilateral du coefficient 
# de correlation de Pearson entre deux vecteurs. La statistique du test est t.
#
# Parameters of the function:
#    x, y: the two data vectors
#    nperm = number of permutations
#
# Example: test the correlation between two vectors of random numbers N(0,1)
#   x <- rnorm(50,0,1)
#   y <- rnorm(50,0,1)
#   cor.out = corPerm2(x, y, 999)
#   cor.out
# Compare the results to:  cor.test(x,y)
#
# Pedagogic objectives --
# - illustrate the structure of a function in the R language
# - show how to program a "for" loop (command: "for(i in 1:nperm)")
# - show how to write out the results using "cat"
# - show how to program a permutation test
#
#                          Pierre Legendre, October 2005
{
x <- as.matrix(x)
y <- as.matrix(y)
n <- nrow(x)

# The object "temp" produced by "cor.test" contains 9 elements; type
# "summary(temp)" to obtain the list of these elements. Among them,
# r = temp$estimate, t = temp$statistic, d.f. = temp$parameter,
# prob = temp$p.value, limits of the 95% C.I. = temp$conf.int
# By default, cor.test produces a two-tailed parametric test ("two.sided")

temp <- cor.test(x,y)
r.ref <- temp$estimate
t.ref <- temp$statistic

nGT <- 1
for(i in 1:nperm)
   {
   y.perm  <-  sample(y,n)
   temp.perm <- cor.test(x,y.perm)
   t.perm <- temp.perm$statistic
   if( abs(t.perm) >= abs(t.ref) ) nGT <- nGT+1
   }
P <- nGT/(nperm+1)
cat('\nPearson correlation (two-tailed test)','\n')
cat('r =',r.ref,'\n')
cat('t =',t.ref,'\n')
cat('d.f. =',temp$parameter,'\n')
cat('95% C.I. of r = [',temp$conf.int[1],' ',temp$conf.int[2],']','\n')
cat('Prob(param) =',temp$p.value,'\n')
cat('Prob(',nperm,'permutations) =',P,'\n','\n')
return(list(Correlation=r.ref, tStat=t.ref, No.perm=nperm, P.perm=P))
}


corPerm3 <- function(x,y,nperm=999,tail=2)
#
# En - This function computes a permutation test of a Pearson correlation
# coefficient between two data vectors. The test statistic is t.
# One-tailed test in the left-hand tail:  tail = -1
# One-tailed test in the right-hand tail: tail =  1
# Two-tailed test:                        tail =  2
#
# Fr - Cette fonction realise, par permutation, un test du coefficient 
# de correlation de Pearson entre deux vecteurs. La statistique du test est t.
# Test unilateral a gauche: tail = -1  
# Test unilateral a droite: tail =  1  
# Test bilateral:           tail =  2
#
# Parameters of the function:
#    x, y: the two data vectors
#    nperm = number of permutations
#
# Example: test the correlation between two vectors of random numbers N(0,1)
#   x <- rnorm(50,0,1)
#   y <- rnorm(50,0,1)
#   cor.out = corPerm3(x, y, 999, 1)    # one-tailed test in the right-hand tail
#   cor.out
# Compare the results to:  cor.test(x, y, alternative="greater")
#
# Pedagogic objectives --
# - illustrate the structure of a function in the R language
# - show how to program a "for" loop (command: "for(i in 1:nperm)")
# - show how to write out the results using "cat"
# - show how to program a permutation test
#
#                          Pierre Legendre, October 2005
{
x <- as.matrix(x)
y <- as.matrix(y)
n <- nrow(x)

# The object "temp" produced by "cor.test" contains 9 elements; type
# "summary(temp)" to obtain the list of these elements. Among them,
# r = temp$estimate, t = temp$statistic, d.f. = temp$parameter,
# prob = temp$p.value, limits of the 95% C.I. = temp$conf.int
# By default, cor.test produces a two-tailed parametric test ("two.sided")

if((tail != -1) & (tail != 1) & (tail != 2)) 
   {stop ("Incorrect value for parameter 'tail'")}
if(tail == -1) temp <- cor.test(x, y, alternative = "less")
if(tail ==  1) temp <- cor.test(x, y, alternative = "greater")
if(tail ==  2) temp <- cor.test(x, y, alternative = "two.sided")
r.ref <- temp$estimate
t.ref <- temp$statistic

nGT <- 1
for(i in 1:nperm)
   {
   y.perm  <-  sample(y,n)
   temp.perm <- cor.test(x,y.perm)
   t.perm <- temp.perm$statistic
   if(tail == -1) if(t.perm <= t.ref) nGT <- nGT+1
   if(tail ==  1) if(t.perm >= t.ref) nGT <- nGT+1
   if(tail ==  2) if( abs(t.perm) >= abs(t.ref) ) nGT <- nGT+1
   }
P <- nGT/(nperm+1)
cat('\nPearson correlation','\n')
if(tail == -1) cat('one-tailed test, left-hand tail','\n')
if(tail ==  1) cat('one-tailed test, right-hand tail','\n')
if(tail ==  2) cat('two-tailed test','\n')
cat('r =',r.ref,'\n')
cat('t =',t.ref,'\n')
cat('d.f. =',temp$parameter,'\n')
cat('95% C.I. of r = [',temp$conf.int[1],' ',temp$conf.int[2],']','\n')
cat('Prob(param) =',temp$p.value,'\n')
cat('Prob(',nperm,'permutations) =',P,'\n','\n')
return(list(Correlation=r.ref, tStat=t.ref, No.perm=nperm, P.perm=P))
}
