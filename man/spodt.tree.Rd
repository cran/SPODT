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
data(dataMALARIA)
#Example : number of malaria episodes per child at each household,
           #from November to December 2009, Bandiagara, Mali.
#Copyright: Pr Ogobara Doumbo, MRTC, Bamako, Mali. email: okd@icermali.org
sp<-spodt(z~loc+x+y, ~1, data=dataMALARIA, weight=TRUE, graft=0.1,
          level.max=7, min.parent=25, min.child=2, rtwo.min=0.01)
spodt.tree(sp)


}

\keyword{SPODT}
\keyword{Spatial Oblique Decision Tree}
\keyword{spatial partitioning}
