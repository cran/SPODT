\name{spodt.map}
\alias{spodt.map}
\title{Mapping the spatial classification issued from SPODT
}
\description{
This graphical function provides the classified map issued from the \code{spodt} function.
}
\usage{
spodt.map(object, x = NULL, y = NULL, map.col = FALSE, col.sgmts = "black")
}
\arguments{
  \item{object}{an object of class \code{spodt}, usually, a result of a call to \code{spodt}.
}
  \item{x, y}{the initial coordinates of each point location.
}
  \item{map.col}{a logical argument coloring point locations if TRUE.
}
  \item{col.sgmts}{an string argument coloring the boudaries between spatial classes.
}
}

\value{The classified map issued from the \code{spodt} function.
}
\references{
\itemize{
\item{Gaudart J, Poudiougou B, Ranque S, Doumbo O. Oblique decision trees for spatial pattern detection: optimal algorithm and application to malaria risk. BMC Medical Research Methodology 2005;5:22}
\item{Gaudart J, Giorgi R, Poudiougou B, Toure O, Ranque S, Doumbo O, Demongeot J. Detection de clusters spatiaux sans point source predefini: utilisation de cinq methodes et comparaison de leurs resultats. Revue d'Epidemiologie et de Sante Publique 2007;55(4):297-306}
\item{Fichet B, Gaudart J, Giusiano B. Bivariate CART with oblique regression trees. International conference of Data Science and Classification, International Federation of Classification Societies, Ljubljana, Slovenia, July 2006.} 
}
}
\author{
Jean Gaudart, Nathalie Graffeo, Guillaume Barbet, Bernard Fichet, Roch Giorgi (Aix-Marseille University)
}


\seealso{
\code{\link{spodt}}, \code{\link{spodt.tree}}, \code{\link{test.spodt}} 
}

\examples{
data(dataExample)
#Example : split the area without covariable analysis
sp<-spodt(dataExample,weight=TRUE,graft=0.05)

spodt.map(sp,dataExample$x,dataExample$y,map.col=TRUE)


}

\keyword{SPODT}
\keyword{Spatial Oblique Decision Tree}
\keyword{spatial partitioning}
