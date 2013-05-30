\name{test.spodt}
\alias{test.spodt}
\title{Monte Carlo hypothesis test of the SPODT classification

}
\description{The \code{test.spodt} function provides Monte Carlo hypothesis test of the final classification issued from the \code{spodt} function. This function performs simulations of the specified null hypothesis and the classication of each simulated data set, using the same rules than the observed dataset classification.

}
\usage{
test.spodt(R2.obs, data, rdist, par.rdist, nb.sim, qt.fact = NULL, ql.fact = NULL,
           weight, graft, level.max, min.parent, min.child, rtwo.min)
}

\arguments{
  \item{R2.obs}{the \code{R2global} issued from the previous \code{spodt} final classification of the observed dataset. Specified as a numerical value between 0 and 1.
}
  \item{data}{a data frame containing the variable and presented in a specified format: first column for location id (loc), second and third column for geographic coordinates (x and y), 4th column for the dependant variable (z), and the following columns for covariates including time
}
  \item{rdist}{a description of the distribution of the dependent variable under the null hypothesis. This can be a character string naming a random generation of a specified distribution, such as \code{"rnorm"}(Gaussian distribution), \code{"rpois"} (Poisson distribution), \code{"rbinom"} (binomial distribution), \code{"runif"} (uniform distribution) ... .
}
  \item{par.rdist}{a list of the parameters needed for the random generation, depending on the null hypothesis distribution, such as \code{c(n,mean,sd)} (Gaussian distribution), \code{c(n,lambda)} (Poisson distribution), \code{c(n,size,prob)} (binomial distribution), \code{c(n,min,max)} (uniform distribution) ... .
}
  \item{nb.sim}{the number of simulation, specified as a positive integer.
}
  \item{ql.fact}{
an optional list of quantitative variables to be split according to a non oblique algorithm. Should be \code{NULL} or specified as a character string. Time factor can be listed as a quantitative factor.
}
  \item{qt.fact}{
an optional list of qualitative variables to be split according to a non oblique algorithm. Should be \code{NULL} or specified as a character string.
}
  \item{weight}{
logical value indicating whether the interclass variances should be weighted or not.
}
  \item{graft}{
if not equals to FALSE, a numerical value between 0 and 1 indicating the minimum improvement of \code{R2global} requires to grafted the final classes.
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

\seealso{\code{\link{spodt}}, \code{\link{spodt.tree}}, \code{\link{spodt.map}} 
}
\examples{
#Example 1:
data(dataExample)
summary(dataExample)
sp<-spodt(dataExample,weight=TRUE,graft=0.05, level.max=5, min.parent=10,
    min.child=5, rtwo.min=0.001)

#to test the previous split using Monte-Carlo approach, and hypothesing a
    #uniform distribution of the dependant variable through the area
test.spodt(sp@R2, dataExample, "runif", c(nrow(dataExample),0,10), 5, weight=TRUE, graft=0.05,
           level.max=5, min.parent=10, min.child=5, rtwo.min=0.001)
#
#Example 2:
#to simulate a dataset with 600 locations
tnbObs <- 600
tloc <- c(rep((10+tnbObs/3):11,2),(20+tnbObs/3+1):(20+2*tnbObs/3))
tx <- runif(tnbObs/3, -1, 1)
ty <- runif(tnbObs/3, -1, 1)
tx <- c(tx,tx,tx)
ty <- c(ty,ty,1/2*ty)
tz <- runif(tnbObs, 0, 1)
tz[which(tx < 0.5 & tx > 0 & ty < 0.5 & ty > 0)] <- 10
tz[which((tx < -0.7  &  ty < -0.7))] <- 10
tdata <- data.frame(tloc,tx,ty,tz)
colnames(tdata)[1:4] <- c("loc", "x", "y", "z")
summary(tdata)
#to split the area:
sp<-spodt(tdata,weight=TRUE,graft=0.05)
#the warning "root is a leaf" tells that no split can be provided by the
    #spodt function according to the splitting parameters
#
#to test the previous split using Monte-Carlo approach, and hypothesing a
    #uniform distribution of the dependant variable through the area
test.spodt(sp@R2, tdata, "runif", c(tnbObs,0,10), 5, weight=TRUE, graft=0.05,
           level.max=5, min.parent=10, min.child=5, rtwo.min=0.001)
}
\keyword{SPODT}
\keyword{Spatial Oblique Decision Tree}
\keyword{spatial partitioning}
