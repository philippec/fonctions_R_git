\name{distgeo}
\alias{geodist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~ geodist - GEOgraphic DISTances ~~ }
\description{
  ~~ 	geodist computes Great Circle Distances (along the Earth's curvature) between
  		a series of of latitude/longitude pairs of points. The resulting distance
  		matrix can be expressed in Radians (rad), Degrees (deg), Nautical knots (kt),
  		Miles (mi) or Kilometers (km). ~~
}
\usage{
distgeo(latitude, longitude, conversion = c("rad", "deg", "kt", "mi", "km"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{latitude}{ ~~ Vector containing latitude of points (in decimal degrees) ~~ }
  \item{longitude}{ ~~ Vector containing longitude of points (in decimal degrees) ~~ }
  \item{conversion}{ ~~ Unit for the resulting distance matrix ~~ }
}
\details{
  ~~ There should be the same number of points in latitude and longitude vectors ~~
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
\references{ ~put references to the literature/web site here ~ }
\author{ ~~ Alain Vaudor, Philippe Casgrain and Kalle Tuomisto ~~ }
\note{ ~~ Original algorithm by A. Vaudor, with modifications by K. Tuomisto. ~~ 

 ~Make other sections like Warning with \section{Warning }{....} ~
}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(latitude, longitude, conversion = c("rad", "deg", "kt", "mi", "km"))
{
	# Arret si latitude et longitude sont de longueur differentes
	if (length(latitude) != length(longitude))
	{
		stop("latitude et longitude bizarre")
	}
	
	if (!is.numeric(latitude))
	{
		stop("latitude n'est pas numerique")
	}

	if (!is.numeric(longitude))
	{
		stop("longitude n'est pas numerique")
	}

	# Number of coordinates
	ncoords <- length(latitude)
	
	# Matrice de resultats
	res.mat <- matrix(nrow=ncoords, ncol=ncoords)
	
	# Calcul de quelques constantes
	rad <- 360 / (2 * pi)
	degconv <- rad / 60
	noeudconv <- 60 * rad
	milleconv <- rad * 1.1516 * 60
	kmconv <- 1.853 * 60 * rad
	epsilon <- 10^-10
	
	latitude <- latitude / rad
	longitude <- longitude / rad
	
	# deux boucles pour remplir la matrice res.mat
	for (i in 1:ncoords)
	{
		for (j in 1:ncoords)
		{
			# Definir les morceaux
			lat1 <- latitude[i]
			lat2 <- latitude[j]
			long1 <- longitude[i]
			long2 <- longitude[j]

			# Calcul des differences de longitudes entre tous les sites
			long1 <- long2 - long1

			# Calcul de latitude entre tous les sites
			lat1 <- (sin(lat1) * sin(lat2)) + (cos(lat1) * cos(lat2) * cos(long1))
			
			# Tous les lat1 se trouvant entre -10^-10 et 10^-10 sont consideres comme des 0
			if (lat2 < epsilon & lat2 > -epsilon)
			{
				lat2<-0
			}
			
			# Si lat1 est plus grand que 0
			if (abs(lat1) > 0)
			{
				lat2 <- 1 - lat1 * lat1
				
				# Si lat2 est plus grand que 0
				if(lat2 > 0)
				{
					lat2 <- atan(sqrt(lat2) / lat1)
				}
				else
				{
					lat2 <- 0
				}
				
				# Si lat1 est negatif
				if (lat1 < 0) 
				{
					lat2 <- pi + lat2
				}
			
				# Si on fait la conversion en radians...
				if(conversion == "rad")
				{
					res.mat[i,j] <- lat2
				}

				# Si on fait la conversion en degres...
				if (conversion == "deg")
				{
					res.mat[i,j] <- lat2 * degconv
				}

				# Si on fait la conversion en noeuds...
				if (conversion == "kt")
				{
					res.mat[i,j] <- lat2 * noeudconv
				}

				# Si on fait la conversion en milles...
				if (conversion == "mi")
				{
					res.mat[i,j] <- lat2 * milleconv
				}

				# Si on fait la conversion en km...
				if (conversion == "km")
				{
					res.mat[i,j] <- lat2 * kmconv
				}
			}
			else
			{
				# Si on fait la conversion en radians...
				if (conversion == "rad")
				{
					res.mat[i,j] <- pi/2
				}

				# Si on fait la conversion en degres...
				if(conversion == "deg")
				{
					res.mat[i,j] <- 90
				}

				# Si on fait la conversion en noeuds...
				if(conversion == "kt")
				{
					res.mat[i,j] <- 5400
				}

				# Si on fait la conversion en milles...
				if (conversion == "mi")
				{
					res.mat[i,j] <- 6218.64
				}

				# Si on fait la conversion en km...
				if (conversion == "km")
				{
					res.mat[i,j] <- 10006.2
				}
			}
		}
	}

	# Transformer res.mat en objet de classe dist
	res <- as.dist(res.mat)
	
	# Renvoyer l'objet res
	return(res)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
