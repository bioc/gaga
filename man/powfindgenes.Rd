\name{powfindgenes}
\alias{powfindgenes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Power computations for differential expression }
\description{
  \code{powfindgenes} evaluates the posterior expected number of true positives
  (e.g. true gene discoveries) if one were to obtain an additional batch
  of data. It uses a gaga model fit based on a set of pilot data.
}
\usage{
powfindgenes(gg.fit, x, groups, batchSize = 1, fdrmax = 0.05, genelimit, v0thre = 1, B = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gg.fit}{GaGa or MiGaGa fit (object of type \code{gagafit}, as returned by \code{fitGG}). }
  \item{x}{\code{ExpressionSet}, \code{exprSet}, data frame or matrix
    containing the gene expression measurements used to fit the model.}
  \item{groups}{If \code{x} is of type \code{ExpressionSet} or
    \code{exprSet}, \code{groups} should be the name of the column
    in \code{pData(x)} with the groups that one wishes to compare. If
    \code{x} is a matrix or a data frame, \code{groups} should be a
    vector indicating to which group each column in x
    corresponds to.}
  \item{batchSize}{ Number of additional samples to obtain per group. }
  \item{fdrmax}{Upper bound on FDR.}.
  \item{genelimit}{ Only the \code{genelimit} genes with the lowest
    probability of being equally expressed across all groups will be
    simulated. Setting this limit can significantly increase the
    computational speed. }
  \item{v0thre}{ Only genes with posterior probability of being equally
    expressed < \code{v0thre} will be simulated.  Setting this limit can significantly increase the
    computational speed.}
  \item{B}{ Number of simulations from the GaGa predictive distribution
    to be used to estimate the posterior expected number of true positives. }
}
\details{
  The routine simulates data from the posterior predictive distribution
  of a GaGa model. That is, first it simulates parameter values (differential
  expression status, mean expression levels etc.) from the posterior
  distribution. Then it simulates data using Gamma distributions and the
  parameter values drawn from the posterior.
  Finally the simulated data is used to determine the differential status
  of each gene, controlling the Bayesian FDR at the \code{fdrmax} level,
  as implemented in \code{findgenes}.
  As the differential expression status is known for each gene, one can
  evaluate the number of true discoveries in the reported gene list.
}
\value{
  \item{m}{Posterior expected number of true positives (as estimated by
    the sample mean of \code{B} simulations)}
  \item{s}{Standard error of the estimate i.e. SD of the simulations/sqrt(B)}
}
\references{ Rossell D. GaGa: a simple and flexible hierarchical model for microarray
data analysis. \url{http://rosselldavid.googlepages.com}. }
\author{ David Rossell }
\seealso{ \code{\link{findgenes}}, \code{\link{fitGG}},
  \code{\link{parest}}. See \code{\link{powclasspred}} for
power calculations for sample classification. }
\examples{
#Simulate data and fit GaGa model
set.seed(1)
x <- simGG(n=20,m=2,p.de=.5,a0=3,nu=.5,balpha=.5,nualpha=25)
gg1 <- fitGG(x,groups=1:2,method='EM')
gg1 <- parest(gg1,x=x,groups=1:2)

#Expected nb of TP for 1 more sample per group
powfindgenes(gg1,x=x,groups=1:2,batchSize=1,fdrmax=.05)$m

#Expected nb of TP for 10 more samples per group
powfindgenes(gg1,x=x,groups=1:2,batchSize=10,fdrmax=.05)$m
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }
\keyword{ models }
