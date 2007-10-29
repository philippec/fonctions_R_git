broken.stick <- function(p)
#
# Compute the expected values of the broken-stick distribution for 'p' pieces.
#
# Example: broken.stick.out.20 = broken.stick(20)
#
#             Pierre Legendre, April 2007
{
result = matrix(0,p,2)
colnames(result) = c("j","E(j)")
for(j in 1:p) {
   E = 0
   for(x in j:p) E = E+(1/x)
   result[j,1] = j
   result[j,2] = E/p
   }
return(result)
}
