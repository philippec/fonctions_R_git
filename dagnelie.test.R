dagnelie.test <- function(X, graph=FALSE)
#
# Performs Dagnelie's test of multinormality on a table X of n objects (rows) 
# and p variables (columns).
#
# Arguments of the function --
#
# X : Multivariate data table (object class: matrix or data.frame).
# graph = TRUE  : Plot a histogram of the Mahalanobis distances between the
#                 objects and the multidimensional mean of all objects.
#       = FALSE : Do not plot a histogram of the Mahalanobis distances
#
# Details --
#
# Dagnelie’s test of multivariate normality based on the Shapiro-Wilk test of 
# normality of Mahalanobis generalized distances is invalid for univariate data 
# (type I error rate too high). Numerical simulations by D. Borcard showed that 
# the test had correct levels of type I error for values of n between 3p and 
# 7.5p, where p is the number of variables in the data table (simulations with 1 
# ≤ p ≤ 50). Outside that range of n values, the results were too liberal, 
# meaning that the test rejected too often the null hypothesis of normality. For 
# p = 2, the simulations showed the test to be valid for 6 ≤ n ≤ 11. If H0 is 
# not rejected in a situation where the test is too liberal, the result is 
# trustworthy.
#
# References --
#
# Dagnelie, P. 1975. L'analyse statistique a plusieurs variables. 
#    Les Presses agronomiques de Gembloux, Gembloux, Belgium.
#
# Legendre, P. and L. Legendre. 1998. Numerical ecology, 2nd English
#    edition. Elsevier, Amsterdam, The Netherlands.
#
# Examples --
#
# Example 1: 2 variables, n = 100
# mat2.dag <- matrix(rnorm(200),100,2)
# (dag2.out <- dagnelie.test(mat2.dag, graph=TRUE))
#
# Example 2: 10 variables, n = 50
# mat10.dag <- matrix(rnorm(500),50,10)
# (dag10.50.out <- dagnelie.test(mat10.dag, graph=TRUE))
#
# Example 3: 10 variables, n = 100
# mat10.dag <- matrix(rnorm(1000),100,10)
# (dag10.100.out <- dagnelie.test(mat10.dag, graph=TRUE))
#
# License: GPL-2 
# Authors: Daniel Borcard and Pierre Legendre, 2007, 2011
{
X <- as.matrix(X)
n <- nrow(X)
p <- ncol(X)

# Compute multidimensional mean vector of all objects
xbar <- apply(X,2,mean)

# Compute inverse of the dispersion matrix
invS <- solve(cov(X))

# Mahalanobis distances between the objects and the multidimensional mean vector
# of all objects (Legendre & Legendre 1998, eq. 4.54 p.184)
D <- as.vector(rep(0,n))
for(i in 1:n) {
   temp <- as.matrix(X[i,]-xbar)
   D[i] <- sqrt(t(temp) %*% invS %*% temp)
   }

if(graph == TRUE) hist(D)

# Shapiro-Wilk test on D
multinorm <- shapiro.test(D)

# Warning messages
note2 <- "The result is trustworthy if H0 is not rejected."
if(p == 1) {
	note <- "Test too liberal for univariate data."
	} else if(p > 50) {
	note = NA 
	note2 = NA 
	} else {
	if((n >= 3*p) & (n <= 7.5*p)) {
		note = "Test result trustworthy, n is between 3*p and 7.5*p"
		note2 <- NA
		} else if(n < 3*p) {
		note = "Test too liberal, n < 3*p"
		} else {
		note = "Test too liberal, n > 7.5*p"
		}
	if(p==2) {
		if(n < 6) note = "Test too liberal, p = 2, n < 6"
		if(n > 11) note = "Test too liberal, p = 2, n > 11"
		}
	}

return(list(Shapiro.Wilk=multinorm, n=n, p=p, Note=note, Note2=note2))
}
