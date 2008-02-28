\name{checkfit}
\alias{checkfit}
\alias{checkfit.gagafit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Check goodness-of-fit of GaGa and MiGaGa models}
\description{
  Produces plots to check fit of GaGa and MiGaGa model. Compares
  observed data with posterior predictive distribution of the model.
  Can also compare posterior distribution of parameters with method of moments estimates.
}
\usage{
checkfit(gg.fit, x, groups, type='data', \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gg.fit}{GaGa or MiGaGa fit (object of type \code{gagafit}, as
    returned by \code{fit.gg}).}
  \item{x}{\code{ExpressionSet}, data frame or matrix
    containing the gene expression measurements used to fit the model.}
  \item{groups}{If \code{x} is of type \code{ExpressionSet},
    \code{groups} should be the name of the column in \code{pData(x)}
    with the groups that one wishes to compare. If \code{x} is a matrix
    or a data frame, \code{groups} should be a vector indicating to which
    group each column in x corresponds to.}
  \item{type}{\code{data} checks marginal density of the data;
    \code{shape} checks shape parameter; \code{mean} checks mean
    parameter; \code{shapemean} checks the joint of shape and mean parameters}
  \item{\dots}{Other arguments to be passed to \code{plot}.}
}
\details{
  The routine generates random draws from the posterior and posterior
  predictive distributions, fixing the hyper-parameters at their
  estimated value (posterior mean if model was fit with
  \code{method=='Bayes'} or maximum likelihood estimate is model was fit
  with \code{method=='EBayes'}).
}
\value{
 Produces a plot.
}
\references{Rossell D. GaGa: a simple and 
flexible hierarchical model for microarray
data analysis. \url{http://rosselldavid.googlepages.com}.}
\author{David Rossell}
\note{Posterior and posterior predictive checks can lack sensitivity to
  detect model misfit, since they are susceptible to over-fitting. An
  alternative is to perform prior predictive checks by generating
  parameters and data with \code{sim.gg}.
}
\seealso{\code{\link{sim.gg}} to simulate samples from the
  prior-predictive distribution, \code{\link{simnewsamples}} to generate parameters and
  observations from the posterior predictive, which is useful to check
  goodness-of-fit individually a desired gene.}
\examples{

}
\keyword{distribution}
\keyword{models}
