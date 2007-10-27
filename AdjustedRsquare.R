AdjustedRsquare <- function(Rsquare, n, m, origin=1)
#
#   This function computes Ezekiel's (1930) adjusted R-square for multiple 
#   regression or the adjusted bimultivariate redundancy statistic for RDA,
#   using as input:
#   
#      Rsquare = R-square obtained by regression or RDA
#      n = number of observations
#      m = number of explanatory variables
#
#   Two options for the adjusted R-square are offered: 
#      origin=0 : in (multiple) regression or RDA through the origin 
#      origin=1 : in ordinary multiple regression or RDA
#
# Example of use: RsquareAdj(0.395165, 125, 9)
#
#           © Pierre Legendre, July 2005
{ 
  if(origin == 0) {    # Regression through the origin
     DFtotal = n
     DFresidual = n-m
     }
  if(origin == 1) {    # Ordinary multiple regression
     DFtotal = n-1
     DFresidual = n-1-m
     }
  adjusted = 1-((1-Rsquare)*DFtotal/DFresidual)
  cat('Rsquare =',Rsquare,'  n =',n,'  m =',m,'\n','Adjusted Rsquare =',adjusted,'\n')
# return(adjusted)   ### Use when this function is called by another function
}
