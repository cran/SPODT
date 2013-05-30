\name{spodt}
\alias{spodt}
\title{
Spatial Oblique Decision Tree main function
}
\description{
The \code{spodt} function provides a spatial partionning.
}
\usage{
spodt(data, ql.fact = NULL, qt.fact = NULL, weight = FALSE, graft = FALSE,
      level.max = 5, min.parent = 10, min.child = 5, rtwo.min = 0.001)
}

\arguments{
  \item{data}{
a data frame containing the variables and presented in a specified format: first column for location ID (loc), second and third column for geographic coordinates (x and y), 4th column for the dependant variable (z), and the following columns for covariates (V1 to Vk) including time (t)
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
The spodt function computes an object of class spodt with the different components of the classification tree, i.e. i) at each step: the point locations within each class, the \code{R2} coefficients of the splitting line; ii) global results: the \code{R2global} (\code{object@R2}), the final partition (\code{object@partition}) including the graft results.
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
\code{\link{spodt.tree}}, \code{\link{spodt.map}}, \code{\link{test.spodt}} 
}
\examples{
data(dataExample)
summary(dataExample)
#Example 1: split the area without covariable analysis
sp<-spodt(dataExample,weight=TRUE,graft=0.05)
sp
sp@R2
#Example 2: split the area with the quantitative covariate
sp<-spodt(dataExample, weight=FALSE, ql.fact = NULL, qt.fact = "qt.fact",
    graft=FALSE, level.max=7, min.parent=5, min.child=3, rtwo.min=0.0001)
sp
sp@R2

#the warning "root is a leaf" tells that no split can be provided by the
    # spodt function according to the splitting parameters
}

\keyword{SPODT}
\keyword{Spatial Oblique Decision Tree}
\keyword{spatial partitioning}
