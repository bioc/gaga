\name{powclasspred}
\alias{powclasspred}
\alias{powclasspred.gagafit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Expected probability that a future sample is correctly classified. }
\description{
  Estimates posterior expected probability that a future sample is correctly
  classified when performing class prediction. The estimate is obtained
  via Monte Carlo simulation from the posterior predictive.
}
\usage{
powclasspred(gg.fit, x, groups, prgroups, v0thre=1, ngene=100, B=100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gg.fit}{GaGa or MiGaGa fit (object of type \code{gagafit}, as returned by \code{fit.gg}). }
  \item{x}{\code{ExpressionSet}, \code{exprSet}, data frame or matrix
    containing the gene expression measurements used to fit the model.}
  \item{groups}{If \code{x} is of type \code{ExpressionSet} or
    \code{exprSet}, \code{groups} should be the name of the column
    in \code{pData(x)} with the groups that one wishes to compare. If
    \code{x} is a matrix or a data frame, \code{groups} should be a
    vector indicating to which group each column in x
    corresponds to.}
  \item{prgroups}{Vector specifying prior probabilities for each
    group. Defaults to equally probable groups.}
  \item{v0thre}{Only genes with posterior probability of being equally
    expressed below \code{v0thre} are used.}
  \item{ngene}{Number of genes to use to build the classifier. Genes
    with smaller probability of being equally expressed are selected
    first.}
  \item{B}{Number of Monte Carlo samples to be used.}
}
\details{
  The routine simulates future samples (microarrays) from the posterior
  predictive distribution of a given group (e.g. control/cancer).
  Then it computes the posterior probability
  that the new sample belongs to each of the groups
  and classifies the sample into the group with
  highest probability. This process is repeated \code{B} times, and the
  proportion of correctly classified samples is reported for each
  group. The standard
  error is obtained via the usual normal approximation (i.e. SD/B).
  The overall probability of correct classification is also provided
  (i.e. for all groups together), but using a more efficient variant of
  the algorithm. Instead of reporting the observed proportion of
  correctly classified samples, it reports the expected proportion of
  correctly classified samples (i.e. the average posterior probability
  of the class that the sample is assigned to).
}
\value{
List with components:
  \item{ccall}{Estimated expected probability of correctly classifying a
    future sample.}
  \item{seccall}{Estimated standard error of \code{ccall}.}
  \item{ccgroup}{Vector with the estimated probability of correctly
    classifying a sample from each group.}
  \item{segroup}{Estimated standard error of \code{ccgroup}.}
}
\references{ Rossell D. GaGa: a simple and 
flexible hierarchical model for microarray
data analysis. \url{http://rosselldavid.googlepages.com}. }
\author{ David Rossell }
\seealso{ \code{\link{classpred}}, \code{\link{fit.gg}},
  \code{\link{parest}} }
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }
\keyword{ models }
