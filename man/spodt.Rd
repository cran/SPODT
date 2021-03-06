\name{spodt}
\alias{spodt}
\title{
Spatial Oblique Decision Tree main function
}
\description{
The \code{spodt} function provides a spatial partionning.
}
\usage{
spodt(formula, data, weight = FALSE, graft = 0,
      level.max = 5, min.parent = 10, min.child = 5, rtwo.min = 0.001)
}

\arguments{
  \item{formula}{
a \link{formula}, with a response but no interaction terms. The left hand side has to contain the quantitative response variable. The right hand side should contain the quantitative and qualitative variables to be split according to a non oblique algorithm. For single spatial analysis (with no cofactor) the right hand side should be ~1.
}
  \item{data}{
a \code{SpatialPointsDataFrame} containing the coordinates and the variables. \code{spodt} needs planar coordinates. Geographic coordinates have to be projected. Otherwise, euclidian coordinates can be used.
}
  \item{weight}{
logical value indicating whether the interclass variances should be weighted or not.
}
  \item{graft}{
if not equals to 0, a numerical value in ]0;1] indicating the minimal modification of \code{R2global} requires to grafted the final classes.
}
  \item{level.max}{
the maximal level of the regression tree above which the splitting algorithm is stopped.
}
  \item{min.parent}{
the minimal size of a node below which the splitting algorithm is stopped.
}
  \item{min.child}{
the minimal size of the children classes below which the split is refused and algorithm is stopped.
}
  \item{rtwo.min}{
the minimal value of \code{R2} above which the node split is refused and algorithm is stopped. Specified as a numerical value between 0 and 1.
}
}

\value{
The spodt function computes an object of class \code{spodt} with the different components of the classification tree, i.e. i) at each step: the point locations within each class, the \code{R2} coefficients of the splitting line; ii) global results: the \code{R2global} (\code{object@R2}), the final partition (\code{object@partition}) including the graft results.
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
\code{\link{spodt.tree}}, \code{\link{spodtSpatialLines}}, \code{\link{test.spodt}} 
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
sp
sp@R2

#the warning "root is a leaf" tells that no split can be provided by the
    # spodt function according to the splitting parameters
}

\keyword{SPODT}
\keyword{Spatial Oblique Decision Tree}
\keyword{spatial partitioning}
