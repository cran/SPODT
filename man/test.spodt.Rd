\name{test.spodt}
\alias{test.spodt}
\title{Monte Carlo hypothesis test of the SPODT classification

}
\description{The \code{test.spodt} function provides Monte Carlo hypothesis test of the final classification issued from the \code{spodt} function. This function performs simulations of the specified null hypothesis and the classification of each simulated data set, using the same rules than the observed dataset classification.

}
\usage{
test.spodt(formula, data, R2.obs,  rdist, par.rdist, nb.sim, 
                   weight=FALSE, graft=0, level.max=5, min.parent=10, 
				   min.child=5, rtwo.min=0.001)
}

\arguments{
  \item{formula}{
a \link{formula}, with a response but no interaction terms. The left hand side has to contain the quantitative response variable. The right hand side should contain the quantitative and qualitative variables to be split according to a non oblique algorithm. For single spatial analysis (with no cofactor) the right hand side should be ~1.
}
  \item{data}{
a \code{SpatialPointsDataFrame} containing the coordinates and the variables. \code{spodt} needs planar coordinates. Geographic coordinates have to be projected. Otherwise, euclidian coordinates can be used.
}
  \item{R2.obs}{the \code{R2global} issued from the previous \code{spodt} final classification of the observed dataset. Specified as a numerical value between 0 and 1.
}
  \item{rdist}{a description of the distribution of the dependent variable under the null hypothesis. This can be a character string naming a random generation of a specified distribution, such as \code{"rnorm"}(Gaussian distribution), \code{"rpois"} (Poisson distribution), \code{"rbinom"} (binomial distribution), \code{"runif"} (uniform distribution) ... .
}
  \item{par.rdist}{a list of the parameters needed for the random generation, depending on the null hypothesis distribution, such as \code{c(n,mean,sd)} (Gaussian distribution), \code{c(n,lambda)} (Poisson distribution), \code{c(n,size,prob)} (binomial distribution), \code{c(n,min,max)} (uniform distribution) ... .
}
  \item{nb.sim}{the number of simulation, specified as a positive integer.
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
The \code{test.spodt} function computes classification trees for simulated dataset. It provides the \code{R2global} empirical distribution under the null hypothesis, compared to the observed \code{R2global}, and a p-value.
}
\references{
\itemize{
\item{Gaudart J, Poudiougou B, Ranque S, Doumbo O. Oblique decision trees for spatial pattern detection: optimal algorithm and application to malaria risk. BMC Medical Research Methodology 2005;5:22}
\item{Gaudart J, Giorgi R, Poudiougou B, Toure O, Ranque S, Doumbo O, Demongeot J. Detection de clusters spatiaux sans point source predefini: utilisation de cinq methodes et comparaison de leurs resultats. Revue d'Epidemiologie et de Sante Publique 2007;55(4):297-306}
\item{Fichet B, Gaudart J, Giusiano B. Bivariate CART with oblique regression trees. International conference of Data Science and Classification, International Federation of Classification Societies, Ljubljana, Slovenia, July 2006.} 
}
}
\author{Jean Gaudart, Nathalie Graffeo, Guillaume Barbet, Bernard Fichet, Roch Giorgi (Aix-Marseille University)
}

\seealso{\code{\link{spodt}}, \code{\link{spodt.tree}}, \code{\link{spodtSpatialLines}} 
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

#to test the previous split using Monte-Carlo approach, and hypothesing a
    #Poisson distribution of the dependant variable through the area
test.spodt(dataMALARIA@data[,2]~1, dataMALARIA, sp@R2, "rpois",
           c(length(dataMALARIA@data$loc),mean(dataMALARIA@data$z)), 10,
		   weight=TRUE, graft=gr, level.max=lmx, min.parent=parm,
		   min.child=childm,rtwo.min=rtw)

#the warning "root is a leaf" tells that no split can be provided by the
    #spodt function according to the splitting parameters

}
\keyword{SPODT}
\keyword{Spatial Oblique Decision Tree}
\keyword{spatial partitioning}
