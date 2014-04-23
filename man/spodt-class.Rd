\name{spodt-class}
\Rdversion{1.1}
\docType{class}
\alias{spodt-class}

\title{Class \code{"spodt"}
Spatial Oblique Decision Tree object}
\description{
This class of objects is returned by the \code{"spodt"} function to represent a spatial oblique decision tree model.
}

\section{Slots}{
  \describe{
    \item{\code{racine}:}{Object of class \code{"virtual.spodt"} the root of the tree }
    \item{\code{R2}:}{Object of class \code{"numeric"} the explained variance coefficient issued from a split }
    \item{\code{partition}:}{Object of class \code{"vector"} the final partition issued from the \code{"spodt"} function }
    \item{\code{adj}:}{Object of class \code{"matrix"} adjacence matrix between spatial classes in a partition. Used for the grafting function }
    \item{\code{cl.grf}:}{Object of class \code{"matrix"} class matrix after grafting }
    \item{\code{sgmts.grf}:}{Object of class \code{"matrix"} border of each class to be grafted }
	\item{\code{brd}:}{Object of class \code{"matrix"} all spatial borders }
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
\author{
Jean Gaudart, Nathalie Graffeo, Guillaume Barbet, Bernard Fichet, Roch Giorgi (Aix-Marseille University)
}


\seealso{
\code{\link{spodt}}, \code{\link{spodt.tree}}, \code{\link{spodtSpatialLines}}, \code{\link{test.spodt}}
}
\examples{
showClass("spodt")
}
\keyword{classes}
\keyword{SPODT}
