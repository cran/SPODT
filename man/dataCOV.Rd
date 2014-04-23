\name{dataCOV}
\alias{dataCOV}

\docType{data}
\title{
Example data file
}
\description{
Simulated data for SPODT functions. To assess the SpODT algorithm to detect 6 situations evoluting through time: rotated square cluster (2 situations), no cluster (2 situations), and V-border shapes (2 situations).
}
\usage{data(dataCOV)}
\format{
  A data frame with 1800 observations on the following 6 variables (300 locations).
  \describe{
    \item{\code{i}}{a numeric vector}
    \item{\code{x}}{a numeric vector}
    \item{\code{y}}{a numeric vector}
    \item{\code{z}}{a numeric vector}
	\item{\code{V1}}{a numeric vector}
	\item{\code{V2}}{a factor with levels ("sq2", "sq1.5", "sq0", "v0", "v1.5", "v2")}
  }
}
\details{
\itemize{
  \item{\code{i}}{: identification of each localization.
}
  \item{\code{x}}{: longitudinal coordinate.
}
  \item{\code{y}}{: latitudinal coordinate.
}
  \item{\code{z}}{: the dependant variable.
}
  \item{\code{V1}}{: a numeric cofactor (1 to 6) to be split according to a non oblique algorithm.
}
  \item{\code{V2}}{:  a qualitative variable to be split according to a non oblique algorithm (6 levels).
}
}
}

\references{
\itemize{
\item{Gaudart J, Graffeo N, Coulibaly D, Barbet G, Rebaudet S, Dessay N, Doumbo O, Giorgi R. SPODT: An R Package to Perform Spatial Partitioning. Journal of Statistical Software 2015;63(16):1-23. \url{http://www.jstatsoft.org/v63/i16/}}
\item{Gaudart J, Poudiougou B, Ranque S, Doumbo O. Oblique decision trees for spatial pattern detection: optimal algorithm and application to malaria risk. BMC Medical Research Methodology 2005;5:22}
\item{Gaudart J, Giorgi R, Poudiougou B, Toure O, Ranque S, Doumbo O, Demongeot J. Detection de clusters spatiaux sans point source predefini: utilisation de cinq methodes et comparaison de leurs resultats. Revue d'Epidemiologie et de Sante Publique 2007;55(4):297-306}
\item{Fichet B, Gaudart J, Giusiano B. Bivariate CART with oblique regression trees. International conference of Data Science and Classification, International Federation of Classification Societies, Ljubljana, Slovenia, July 2006.}
}
}
\examples{
#data(dataCOV)
#dataset<-dataCOV
#coordinates(dataset)<-c("x","y")
#coordinates are planar ones
#Example : split the area with a covariable analysis
#sp<-spodt(dataset@data$z~V1, dataset, weight=FALSE, graft=0, 
#         min.ch=150, min.parent=200, level.max=3, rtwo.min=0.15)
#spodt.tree(sp)

}
\keyword{datasets}
