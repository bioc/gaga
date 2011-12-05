\name{powfindgenes}
\alias{powfindgenes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Power computations for differential expression }
\description{
  \code{powfindgenes} evaluates the posterior expected number of true positives
  (e.g. true gene discoveries) if one were to obtain an additional batch
  of data. It uses either a GaGa or a normal-normal model fit on a pilot
  data set.
}
\usage{
powfindgenes(fit, x, groups, batchSize = 1, fdrmax = 0.05, genelimit,
v0thre = 1, B = 1000, mc.cores=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{GaGa/MiGaGa or normal-normal model fit using pilot data
    \code{x}. It must be an object either of type \code{gagafit} (see
    \code{fitGG}) or \code{nnfit} (see \code{fitNN}).}
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
    to be used to estimate the posterior expected number of true
    positives. }
  \item{mc.cores}{ If \code{multicore} package is available,
    \code{mc.cores} indicates the number of cores to use for parallel
    computing. Currently only used when \code{fit} is of class \code{nnfit}.}
}
\details{
  The routine simulates data from the posterior predictive distribution
  of a GaGa or normal-normal model. That is, first it simulates parameter values (differential
  expression status, mean expression levels etc.) from the posterior
  distribution. Then it simulates data using the
  parameter values drawn from the posterior.
  Finally the simulated data is used to determine the differential status
  of each gene, controlling the Bayesian FDR at the \code{fdrmax} level,
  as implemented in \code{findgenes}.
  As the differential expression status is known for each gene, one can
  evaluate the number of true discoveries in the reported gene list.

  In order to improve speed, hyper-parameters are not re-estimated when
  computing posterior probabilities for the posterior predictive
  simulated data.
}
\value{
  \item{m}{Posterior expected number of true positives (as estimated by
    the sample mean of \code{B} simulations)}
  \item{s}{Standard error of the estimate i.e. SD of the simulations/sqrt(B)}
}
\references{ Rossell D. GaGa: a simple and flexible hierarchical model for microarray
data analysis. \url{http://rosselldavid.googlepages.com}. }
\author{ David Rossell }
\seealso{ \code{\link{findgenes}}, \code{\link{fitGG}}, \code{\link{fitNN}},
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
