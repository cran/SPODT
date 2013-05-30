\name{SPODT-package}
\alias{SPODT-package}
\alias{SPODT}
\docType{package}
\title{
SPatial Oblique Decision Tree
}
\description{
SPODT is a spatial partitioning method based on oblique decision trees, in order to classify study area into zones of different risks, determining their boundaries
}
\details{
\tabular{ll}{
Package: \tab SPODT\cr
Type: \tab Package\cr
Version: \tab 0.6.32\cr
Date: \tab 2013-05-30\cr
License: GPL(>=2)
}
The main functions are spodt(), that provides the spatial classification, spodt.tree(), that provides the partition tree, spodt.map(), that provides the final spatial classification, and test.spodt(), that provides a Monte Carlo test of the final spatial classification.
}
\author{
Jean Gaudart, Nathalie Graffeo, Guillaume Barbet, Bernard Fichet, Roch Giorgi (Aix-Marseille University)

Maintainer: Jean Gaudart <jean.gaudart@univ-amu.fr>

}
\references{
\itemize{
\item{Gaudart J, Poudiougou B, Ranque S, Doumbo O. Oblique decision trees for spatial pattern detection: optimal algorithm and application to malaria risk. BMC Medical Research Methodology 2005;5:22}
\item{Gaudart J, Giorgi R, Poudiougou B, Toure O, Ranque S, Doumbo O, Demongeot J. Detection de clusters spatiaux sans point source predefini: utilisation de cinq methodes et comparaison de leurs resultats. Revue d'Epidemiologie et de Sante Publique 2007;55(4):297-306}
\item{Fichet B, Gaudart J, Giusiano B. Bivariate CART with oblique regression trees. International conference of Data Science and Classification, International Federation of Classification Societies, Ljubljana, Slovenia, July 2006.} 
}

}

\keyword{package}

\examples{

data(dataMALARIA)
#Example : number of malaria episodes per child at each household,
          #from November to December 2009, Bandiagara, Mali.
#Copyright: Pr Ogobara Doumbo, MRTC, Bamako, Mali. email: okd@icermali.org
sp<-spodt(z~loc+x+y, ~1, data=dataMALARIA, weight=TRUE, graft=0.1,
          level.max=7, min.parent=25, min.child=2, rtwo.min=0.01)
spodt.tree(sp)
spodt.map(sp,dataMALARIA$x,dataMALARIA$y,map.col=TRUE)

#the warning "root is a leaf" tells that no split can be provided by the
    #spodt function according to the splitting parameters
}
