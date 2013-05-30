\name{dataBALL0_5}
\alias{dataBALL0_5}

\docType{data}
\title{
Example data file
}
\description{
Simulated data for SPODT functions. To assess the SpODT algorithm to detect a ball shape situation: Clustered data with moderate high values outside a ball shape. 
}
\usage{data(dataBALL0_5)}
\format{
  A data frame with 300 observations on the following 4 variables (300 locations).
  \describe{
    \item{\code{i}}{a numeric vector}
    \item{\code{x}}{a numeric vector}
    \item{\code{y}}{a numeric vector}
    \item{\code{z}}{a numeric vector}
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
}
}

\references{
\itemize{
\item{Gaudart J, Poudiougou B, Ranque S, Doumbo O. Oblique decision trees for spatial pattern detection: optimal algorithm and application to malaria risk. BMC Medical Research Methodology 2005;5:22}
\item{Gaudart J, Giorgi R, Poudiougou B, Toure O, Ranque S, Doumbo O, Demongeot J. Detection de clusters spatiaux sans point source predefini: utilisation de cinq methodes et comparaison de leurs resultats. Revue d'Epidemiologie et de Sante Publique 2007;55(4):297-306}
\item{Fichet B, Gaudart J, Giusiano B. Bivariate CART with oblique regression trees. International conference of Data Science and Classification, International Federation of Classification Societies, Ljubljana, Slovenia, July 2006.}
}
}
\examples{
data(dataBALL0_5)
#Example : split the area without covariable analysis
sp<-spodt(z~i+x+y, ~1, data=dataBALL0_5,weight=TRUE,graft=0.2)
spodt.tree(sp)
spodt.map(sp,dataBALL0_5$x,dataBALL0_5$y,map.col=TRUE)


}
\keyword{datasets}
