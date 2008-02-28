\name{simnewsamples}
\alias{simnewsamples}
\alias{simnewsamples.gagafit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Posterior predictive simulation}
\description{
  Simulates parameters and data from the posterior and posterior
  predictive distributions, respectively, of a GaGa or MiGaGa model.
}
\usage{
simnewsamples(gg.fit, groupsnew, sel, x, groups)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gg.fit}{GaGa or MiGaGa fit (object of type \code{gagafit}, as returned by \code{fit.gg}).}
  \item{groupsnew}{Vector indicating the group that each new sample
    should belong to. \code{length(groupsnew)} is the number of new
    samples that will be generated.}
  \item{sel}{Numeric vector with the indexes of the genes we want to draw new samples for (defaults to all genes). If a logical vector is indicated, it is converted to \code{(1:nrow(x))[sel]}.}
  \item{x}{\code{ExpressionSet}, data frame or matrix
    containing the gene expression measurements used to fit the model.}
  \item{groups}{If \code{x} is of type \code{ExpressionSet},
    \code{groups} should be the name of the column in \code{pData(x)}
    with the groups that one wishes to compare. If \code{x} is a matrix
    or a data frame, \code{groups} should be a vector indicating to which
    group each column in x corresponds to.}
}
\details{
  The shape parameters are actually drawn from a gamma approximation to
  their posterior distribution. The function \code{rcgamma} implements
  this approximation.
}
\value{
  List with the following components:
  \item{xnew}{Matrix of \code{length(nsel)} rows and
    \code{length(groupsnew)} columns with observations drawn from the
    posterior predictive.}
  \item{dnew}{Matrix with same dimensions as \code{xnew} with with expression patterns drawn from the posterior.}
  \item{anew}{Matrix with same dimensions as \code{xnew} with with shape
    parameters drawn from the posterior.}
  \item{lnew}{Matrix with same dimensions as \code{xnew} with with mean
    expression drawn from the posterior.}
}
\references{ Rossell D. GaGa: a simple and 
flexible hierarchical model for microarray
data analysis. \url{http://rosselldavid.googlepages.com}. }
\author{ David Rossell }
\seealso{ \code{\link{checkfit}} for posterior predictive plot,
  \code{\link{sim.gg}} for prior predictive simulation. }
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ distribution }
\keyword{ models }
