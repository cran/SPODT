\name{dataExample}
\alias{dataExample}

\docType{data}
\title{
Example data file
}
\description{
Simulated data for SPODT functions.
}
\usage{data(dataExample)}
\format{
  A data frame with 300 observations on the following 6 variables.
  \describe{
    \item{\code{loc}}{a numeric vector}
    \item{\code{x}}{a numeric vector}
    \item{\code{y}}{a numeric vector}
    \item{\code{z}}{a numeric vector}
    \item{\code{qt.fact}}{a numeric vector}
    \item{\code{ql.fact}}{a factor with levels \code{1} \code{2} \code{3} \code{4} \code{5}}
  }
}
\details{
\itemize{
  \item{\code{loc}}{: identification of each localization.
}
  \item{\code{x}}{: longitudinal coordinate.
}
  \item{\code{y}}{: latitudinal coordinate.
}
  \item{\code{z}}{: the dependant variable.
}
  \item{\code{ql.fact}}{: a quantitative variable to be split according to a non oblique algorithm.
}
  \item{\code{qt.fact}}{: a qualitative variable to be split according to a non oblique algorithm.
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
data(dataExample)
#Example : split the area without covariable analysis
sp<-spodt(dataExample,weight=TRUE,graft=0.05)
spodt.tree(sp)
spodt.map(sp,dataExample$x,dataExample$y,map.col=TRUE)


}
\keyword{datasets}
