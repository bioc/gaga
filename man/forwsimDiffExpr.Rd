\name{forwsimDiffExpr}
\alias{forwsimDiffExpr}
\alias{forwsimDiffExpr.gagafit}
\alias{forwsimDiffExpr.nnfit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Forward simulation for differential expression. }
\description{
  Forward simulation allows to evaluate the expected utility for
  sequential designs. Here the utility is the expected number of true
  discoveries minus a sampling cost.
  The routine simulates future data either from the prior predictive or
  using a set of pilot data and a GaGa or normal-normal
  model fit. At each future time point, it computes a summary statistic
  that will be used to determine when to stop the experiment.
}
\usage{
forwsimDiffExpr(fit, x, groups, ngenes, maxBatch, batchSize, fdrmax = 0.05, genelimit, v0thre = 1, B = 100,
Bsummary = 100, trace = TRUE, randomSeed, mc.cores=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{Either GaGa or MiGaGa fit (object of type \code{gagafit},
    as returned by \code{fitGG}) or Normal-Normal fit (type \code{nnfit}
    returned by \code{fitNN}). }
  \item{x}{\code{ExpressionSet}, \code{exprSet}, data frame or matrix
    containing the gene expression measurements used to fit the model.}
  \item{groups}{If \code{x} is of type \code{ExpressionSet} or
    \code{exprSet}, \code{groups} should be the name of the column
    in \code{pData(x)} with the groups that one wishes to compare. If
    \code{x} is a matrix or a data frame, \code{groups} should be a
    vector indicating to which group each column in x
    corresponds to.}
  \item{ngenes}{Number of genes to simulate data for. If \code{x} is
    specified this argument is set to \code{nrow(x)} and data is
    simulated from the posterior predictive conditional on \code{x}. If
    \code{x} not specified simulation is from the prior predictive. }
  \item{maxBatch}{ Maximum number of batches, i.e. the routine simulates
    \code{batchSize}*\code{maxBatch} samples per group.}
  \item{batchSize}{ Batch size, i.e. number of observations per group to
    simulate at each time point. Defaults to \code{ncol(x)/length(unique(groups))}. }
  \item{fdrmax}{Upper bound on FDR.}
  \item{genelimit}{ Only the \code{genelimit} genes with the lowest
    probability of being equally expressed across all groups will be
    simulated. Setting this limit can significantly increase the
    computational speed. }
  \item{v0thre}{ Only genes with posterior probability of being equally
    expressed < \code{v0thre} will be simulated.  Setting this limit can significantly increase the
    computational speed.}
  \item{B}{ Number of forward simulations. }
  \item{Bsummary}{ Number of simulations for estimating the summary statistic. }
  \item{trace}{For \code{trace==TRUE} iteration progress is displayed.}
  \item{randomSeed}{Integer value used to set random number generator
    seed. Defaults to \code{as.numeric(Sys.time())} modulus 10^6.}
  \item{mc.cores}{ If \code{multicore} package is available,
    \code{mc.cores} indicates the number of cores to use for parallel
    computing. Currently only used when \code{fit} is of class \code{nnfit}.}
}
\details{
To improve computational speed hyper-parameters are not re-estimated as
new data is simulated.
}
\value{
A \code{data.frame} with the following columns:
  \item{simid}{Simulation number.}
  \item{j}{Time (sample size).}
  \item{u}{Expected number of true positives if we were to stop
    experimentation at this time.}
  \item{fdr}{Expected FDR if we were to stop experimentation at this time.}
  \item{fnr}{Expected FNR if we were to stop experimentation at this time.}
  \item{power}{Expected power (as estimated by E(TP)/E(positives)) if we
    were to stop experimentation at this time.}
  \item{summary}{Summary statistic: increase in expected true positives
    if we were to obtain one more data batch.}
}
\references{
  Rossell D., Mueller P. Sequential sample sizes for
  high-throughput hypothesis testing
  experiments. \url{http://sites.google.com/site/rosselldavid/home}.
  
  Rossell D. GaGa: a simple and flexible hierarchical model for microarray
   data analysis. Annals of Applied Statistics, 2009, 3, 1035-1051.
}
\author{ David Rossell. }
\seealso{ \code{\link{fitGG}} for fitting a GaGa model,
  \code{\link{fitNN}} for fitting a normal-normal model,
  \code{\link{seqBoundariesGrid}} for finding the optimal design based
  on the forwards simulation output.
\code{\link{powfindgenes}} for fixed sample size calculations.}
\examples{
#Simulate data and fit GaGa model
set.seed(1)
x <- simGG(n=20,m=2,p.de=.5,a0=3,nu=.5,balpha=.5,nualpha=25)
gg1 <- fitGG(x,groups=1:2,method='EM')
gg1 <- parest(gg1,x=x,groups=1:2)

#Run forward simulation
fs1 <- forwsimDiffExpr(gg1, x=x, groups=1:2,
maxBatch=2,batchSize=1,fdrmax=0.05, B=100, Bsummary=100, randomSeed=1)

#Expected number of true positives for each sample size
tapply(fs1$u,fs1$time,'mean')

#Expected utility for each sample size
samplingCost <- 0.01
tapply(fs1$u,fs1$time,'mean') - samplingCost*(0:2)

#Optimal sequential design
b0seq <- seq(0,20,length=200); b1seq <- seq(0,40,length=200)
bopt <-seqBoundariesGrid(b0=b0seq,b1=b1seq,forwsim=fs1,samplingCost=samplingCost,powmin=0)
bopt <- bopt$opt

plot(fs1$time,fs1$u,xlab='Additional batches',ylab='E(newly discovered DE genes)')
abline(bopt['b0'],bopt['b1'])
text(.2,bopt['b0'],'Continue',pos=3)
text(.2,bopt['b0'],'Stop',pos=1)
}
\keyword{ htest }
\keyword{ design }
