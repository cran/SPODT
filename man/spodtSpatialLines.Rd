\name{spodtSpatialLines}
\alias{spodtSpatialLines}
\title{Boundaries of the spatial classification issued from SPODT
}
\description{
This function creates a SpatialLines object that contains the boundaries of the spatial classification issued from the \code{spodt} function.
}
\usage{
spodtSpatialLines(object, data)
}
\arguments{
  \item{object}{an object of class \code{spodt}, usually, a result of a call to \code{spodt}.
}
  \item{data}{
a \code{SpatialPointsDataFrame} containing the coordinates and the variables. \code{spodt} needs planar coordinates. Geographic coordinates have to be projected. Otherwise, euclidian coordinates can be used.
}
}

\value{A SpatialLines object that contains the boundaries of the spatial classification issued from the \code{spodt} function.
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
\code{\link{spodt}}, \code{\link{spodt.tree}}, \code{\link{test.spodt}} 
}

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

ssp<-spodtSpatialLines(sp,dataMALARIA)
plot(ssp)
points(dataMALARIA,cex=log(dataMALARIA@data$z*10))


}

\keyword{SPODT}
\keyword{Spatial Oblique Decision Tree}
\keyword{spatial partitioning}
