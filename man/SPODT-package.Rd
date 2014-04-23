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
Version: \tab 0.9\cr
Date: \tab 2014-04-23\cr
License: GPL(>=2)
}
The main functions are spodt() that provides the spatial classification, spodt.tree() that provides the partition tree, spodtSpatialLines(), an object of SpatialLines class that contains the final spatial classification, and test.spodt() that provides a Monte Carlo test of the final spatial classification.
}
\author{
Jean Gaudart, Nathalie Graffeo, Guillaume Barbet, Bernard Fichet, Roch Giorgi (Aix-Marseille University)

Maintainer: Jean Gaudart <jean.gaudart@univ-amu.fr>

}
\references{
\itemize{
\item{Gaudart J, Graffeo N, Coulibaly D, Barbet G, Rebaudet S, Dessay N, Doumbo O, Giorgi R. SPODT: An R Package to Perform Spatial Partitioning. Journal of Statistical Software 2015;63(16):1-23. \url{http://www.jstatsoft.org/v63/i16/}}
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
#Copyright: Pr Ogobara Doumbo, MRTC, Bamako, Mali. email: okd[at]icermali.org
coordinates(dataMALARIA)<-c("x","y")
class(dataMALARIA)
proj4string(dataMALARIA)<-"+proj=longlat +datum=WGS84 +ellps=WGS84"
dataMALARIA<-spTransform(dataMALARIA, CRS("+proj=merc +datum=WGS84 +ellps=WGS84"))

gr<-0.07   #graft parameter
rtw<-0.01 #rtwo.min
parm<-25  #min.parent
childm<-2 #min.child
lmx<-7 

sp<-spodt(dataMALARIA@data[,2]~1, dataMALARIA, weight=TRUE, graft=gr, min.ch=childm,
          min.parent=parm, level.max=lmx, rtwo.min=rtw)
		  
spodt.tree(sp)

ssp<-spodtSpatialLines(sp,dataMALARIA)
plot(ssp)
points(dataMALARIA,cex=log(dataMALARIA@data$z*10))

#test.spodt(dataMALARIA@data[,2]~1, dataMALARIA, sp@R2, "rpois",
#           c(length(dataMALARIA@data$loc),mean(dataMALARIA@data$z)), 10,
#		   weight=TRUE, graft=gr, level.max=lmx, min.parent=parm,
#		   min.child=childm,rtwo.min=rtw)

#the warning "root is a leaf" tells that no split can be provided by the
    #spodt function according to the splitting parameters
}
