\name{spodt.tree}
\alias{spodt.tree}
\title{Spatial Classification Tree

}
\description{
This graphical function provides the tree issued from the \code{spodt} function. 
}
\usage{
spodt.tree(object)
}

\arguments{
  \item{object}{
object is an object of class \code{spodt}, usually, a result of a call to spodt.
}
}

\value{
Each step of the classification are presented with main statistics.  For graphical convenience, grafted classes are not presented but only indicated by class id number.
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
\code{\link{spodt}}, \code{\link{spodt.map}}, \code{\link{test.spodt}} 
}
\examples{
data(dataExample)
#Example 1: split the area without covariable analysis
sp<-spodt(dataExample,weight=TRUE,graft=0.05, rtwo.min = 0.0001)
spodt.tree(sp)


#Example 2: split the area with the quantitative covariate
sp<-spodt(dataExample, weight=FALSE, ql.fact = NULL, qt.fact = "qt.fact",
    graft=FALSE, level.max=7, min.parent=5, min.child=3, rtwo.min=0.0001)
spodt.tree(sp)


}

\keyword{SPODT}
\keyword{Spatial Oblique Decision Tree}
\keyword{spatial partitioning}
