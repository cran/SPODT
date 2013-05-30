\name{dataMALARIA}
\alias{dataMALARIA}

\docType{data}
\title{
Example data file
}
\description{
Field observation data. To assess the SpODT algorithm analysing the number of malaria episodes per child at each household, from November to December 2009, Bandiagara, Mali.
Copyright: Pr Ogobara Doumbo, MRTC, Bamako, Mali. email: okd@icermali.org
}
\usage{data(dataMALARIA)}
\format{
  A data frame with 168 observations on the following 4 variables (168 locations).
  \describe{
    \item{\code{loc}}{a numeric vector}
    \item{\code{x}}{a numeric vector}
    \item{\code{y}}{a numeric vector}
    \item{\code{z}}{a numeric vector}
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
data(dataMALARIA)
#Example : number of malaria episodes per child at each household,
          #from November to December 2009, Bandiagara, Mali.
#Copyright: Pr Ogobara Doumbo, MRTC, Bamako, Mali. email: okd@icermali.org
sp<-spodt(z~loc+x+y, ~1, data=dataMALARIA, weight=TRUE, graft=0.1,
          level.max=7, min.parent=25, min.child=2, rtwo.min=0.01)
spodt.tree(sp)
spodt.map(sp,dataMALARIA$x,dataMALARIA$y,map.col=TRUE)


}
\keyword{datasets}
