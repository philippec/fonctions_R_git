\name{forward.sel}
\alias{forward.sel}
\alias{forward.sel.par}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Forward selection with multivariate Y using a parametric model
  or by permutation under reducel model }
\description{
  Performs a forward selection by permutation of residuals under reduced model. Y can be multivariate.
  }
\usage{
forward.sel(Y, X, K = nrow(X) - 1, R2thresh = 0.99, adjR2thresh = 0.99,nperm = 999, 
    R2more = 0.001, alpha = 0.05, Xscale = TRUE, Ycenter = TRUE, Yscale
= FALSE)
forward.sel.par(Y, X, alpha = 0.05, K = nrow(X)-1, R2thresh = 0.99, R2more = 0.001, adjR2thresh = 0.99, Yscale = FALSE, verbose=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{ A matrix of n lines and m columns that contains (numeric) response variables.  }
  \item{X}{ A matrix of n lines and p columns that contains (numeric) explanatory variables.  }
  \item{K}{ This number is the number of variables to be selected in the forward selection. The default setting is one minus the number of row.(See details for more information) }
  \item{R2thresh}{ The number given here is a R2 parameter. If the
    forward selection has a selection of variable which represent the
    number presented in the parameter or higher, after the introduction
    of a variable, the forward selection will stop. The setting of this
    parameter varies from 0.01 to 1. (See details for more information)
  }
  \item{adjR2thresh}{ The number given here is a adjusted R2 parameter. If the forward selection has a selection of variable which represent the number presented in the parameter or higher, after the introduction of a variable, the forward selection will stop. The setting of this parameter varies from 0.01 to 1. (See details for more information) }
  \item{nperm}{ The number of permutation to be done on the forward selection. the default setting is 999 permutation. }
  \item{R2more}{ The number given here is a R2 parameter. If the forward selection gets to a point where the R2 given by a variable is lower than R2more it will stops. The default setting is 0.001. (See details for more information) }
  \item{alpha}{ The number given here is a significance level. If the p-value of a variable is higher than alpha, the procedure stops. The default setting is 0.05. (See details for more information) }
  \item{Xscale}{ This parameter scales the data entered as parameter X. The default setting is TRUE }
  \item{Ycenter}{ This parameter centers the data entered as parameter Y. The default setting is TRUE }
  \item{Yscale}{ This parameter scales the data entered as parameter Y. The default setting is FALSE }
  \item{verbose}{If 'TRUE' more diagnostics are printed. The default setting is TRUE } 
}
\details{
 The forward selection will stop when either K, R2tresh, adjR2tresh,
 alpha and R2more has its parameter reached. The parametric test for the increase in R-square statistic in forward selection, as implemented in the function \code{forward.sel.par}, can be applied as follows.

(a) If Y is univariate, this function implements the standard parametric F-test used in forward selection (FS) in multiple regression.

(b) If Y is multivariate, this function implements FS using the modified
F-test described by Miller and Farr (1971). This test requires that

-- the Y variables be standardized,

  -- the error in the response variables be normally distributed. This condition must be verified by the user.
}
\value{ 
  A dataframe with:
  \item{ variables }{ The names of the variables }
  \item{ order }{ The order of the selection of the variables  }
  \item{ R2 }{ The R2 of the variable selected }
  \item{ R2Cum }{ The cumulative R2 of the variables selected }
  \item{ AdjR2Cum }{ The cumulative adjusted R2 of the variables selected }
  \item{ F }{ The F statistic }
  \item{ pval }{ The P-value statistic }
  }
  \references{ Canoco manual p.49 \cr
Miller, J. K., and S. D. Farr. (1971). Bimultivariate redundancy: a comprehensive measure of interbattery relationship. \emph{Multivariate Behavioral Research}, \bold{6}, 313--324.\cr

  }
\author{ Stephane Dray. For the parametric method, original code of
  Pierre Legendre and Guillaume Blanchet.}
\note{ Not yet implemented for CCA (weighted regression) and with covariables. }

\examples{
    x=matrix(rnorm(30),10,3)
    y=matrix(rnorm(50),10,5)
    
    forward.sel(y,x,nperm=99, alpha = 0.5)
  }
}
\keyword{ multivariate }% at least one, from doc/KEYWORDS
