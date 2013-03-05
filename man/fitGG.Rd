\name{fitGG}
\alias{fitGG}
\alias{fitNN}
\alias{fitNNSingleHyp}
\alias{adjustfitNN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Fit GaGa hierarchical model }
\description{
 fitGG fits GaGa/MiGaGa hierarchical models, either via a fully Bayesian
 approach or via maximum likelihood.

 fitNN fits a normal-normal hierarchical model (wrapper to
 call emfit in package EBarrays with the LNNMV model).
 fitNNSingleHyp is the same as fitNN but only considers the pattern that all
 groups are different from each other.
 
 adjustfitNN corrects a small sample-size bias in the fitNN estimation procedure.
}
\usage{
fitGG(x, groups, patterns, equalcv = TRUE, nclust = 1, method =
 "quickEM", B, priorpar, parini, trace = TRUE)

fitNN(x, groups, patterns, B=20, trace=TRUE)

fitNNSingleHyp(x, groups, B=10, trace=TRUE)

adjustfitNN(fit, pitrue, B=5, nsim=3, mc.cores=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{\code{ExpressionSet}, \code{exprSet}, data frame or matrix
    containing the gene expression measurements used to fit the
    model. For \code{fitNN}, \code{x} should be in log-scale (in
    contrast to the input to \code{emfit}). \code{fitGG} accepts raw and
    log-scale (as long as all log-measurements remain positive). By
    default we recommend using log-scale to reduce the influence of outliers. }
  \item{groups}{If \code{x} is of type \code{ExpressionSet} or
    \code{exprSet}, \code{groups} should be the name of the column
    in \code{pData(x)} with the groups that one wishes to compare. If
    \code{x} is a matrix or a data frame, \code{groups} should be a
    vector indicating to which group each column in x
    corresponds to.}
  \item{patterns}{Matrix indicating which groups are put together under
    each pattern, i.e. the hypotheses to consider for each
    gene. \code{colnames(patterns)} must match the group levels
    specified in \code{groups}.
    Defaults to two hypotheses: null hypothesis of all groups
    being equal and full alternative of all groups being different.
    The function \code{buildPatterns} can be used to construct a matrix
    with all possible patterns.}
  \item{equalcv}{\code{equalcv==TRUE} fits model assuming constant CV across groups. \code{equalcv==FALSE} compares cv as well as mean expression levels between groups}
  \item{nclust}{Number of clusters in the MiGaGa model. \code{nclust}
    corresponds to the GaGa model. }
  \item{method}{ \code{method=='MH'} fits a fully Bayesian model via
    Metropolis-Hastings posterior sampling. \code{method=='Gibbs'} does
    the same using Gibbs sampling. \code{method=='SA'} uses Simulated
    Annealing to find the posterior mode. \code{method=='EM'} finds
    maximum-likelihood estimates via the expectation-maximization
    algorithm, but this is currently only implemented for
    \code{nclust>1}. \code{method=='quickEM'} is a quicker
    implementation that only performs 2 optimization steps (see details).}
  \item{B}{ Number of iterations to fit the model. For \code{method=='MH'} and \code{method=='Gibbs'}, \code{B}
    is the number of MCMC iterations (defaults to 1000). For
    \code{method=='SA'}, \code{B} is the number of iterations in the
    Simulated Annealing scheme (defaults to 200). For
    \code{method=='EM'}, \code{B} is the maximum number of iterations
    (defaults to 20). }
  \item{priorpar}{ List with prior parameter values. It must have
    components \code{a.alpha0,b.alpha0,a.nu,b.nu,a.balpha,b.balpha,a.nualpha,b.nualpha,p.probclus}
    and \code{p.probpat}. If missing they are set to non-informative
    values that are usually reasonable for RMA and GCRMA normalized data.}
  \item{parini}{ list with components \code{a0}, \code{nu},
    \code{balpha}, \code{nualpha}, \code{probclus} and \code{probpat}
    indicating the starting values for the hyper-parameters. If not
    specified, a method of moments estimate is used.}
  \item{trace}{ For \code{trace==TRUE} the progress of the model fitting
    routine is printed.}
  \item{fit}{\code{nnfit} object, as returned by \code{fitNN}}
  \item{pitrue}{Grid of true \code{pi} values for which to evaluate the
    MoM estimation bias. See details.}
  \item{nsim}{Number of datasets to simulate for each \code{pitrue} value}
  \item{mc.cores}{If package \code{multicore} is available,
    \code{mc.cores} specifies the number of cores to be used for
    parallel computing.}
}
\details{
  For GaGa/MiGaGa models, an approximation is used to sample faster from the
  posterior distribution of the gamma shape parameters and to compute
  the normalization constants (needed to evaluate the likelihood). These
  approximations are implemented in \code{rcgamma} and \code{mcgamma}.

  The cooling scheme in \code{method=='SA'} uses a temperature equal to
  \code{1/log(1+i)}, where \code{i} is the iteration number.

  The EM implementation in \code{method=='quickEM'} is a quick EM
  algorithm that usually delivers hyper-parameter estimates very similar
  to those obtained via the slower \code{method=='EM'}. Additionally,
  the GaGa model inference has been seen to be robust to moderate
  changes in the hyper-parameter estimates in most datasets.

  \code{fitNN} is a wrapper to \code{emfit} in package EBarrays with the LNNMV model.
  This procedure estimates hyper-parameters using the method of moments
  (MoM), which typically results in over-estimating the proportion of
  differentially expressed genes, which we denote by pi1. 
  \code{adjustfitNN} corrects this bias by repeatedly simulating from
  the prior predictive of a normal-normal model. Simulations are
  performed for a grid of pi1 values, so that the expected bias can be
  evaluated for each of them. The bias is then modeled as a smooth
  function of pi1 using function \code{gam} from package \code{mgcv}.
  Finally, the value of pi1 is bias-adjusted and the posterior
  probabilities are recomputed using the updated pi1 value.
}
\value{ \code{fitGG} returns an object of class \code{gagafit}, with components
  \item{parest }{Hyper-parameter estimates. Only returned if \code{method=='EBayes'}, for \code{method=='Bayes'} one must call the function \code{parest} after \code{fitGG}}
  \item{mcmc }{Object of class \code{mcmc} with posterior draws for hyper-parameters. Only returned if \code{method=='Bayes'}.}
  \item{lhood}{For \code{method=='Bayes'} it is the log-likelihood evaluated at each MCMC iteration. For \code{method=='EBayes'} it is the log-likelihood evaluated at the maximum.}
  \item{nclust}{Same as input argument.}
  \item{patterns}{Same as input argument, converted to object of class
    \code{gagahyp}.}

  \code{fitNN} returns an analogous object of class \code{nnfit}. The
  component \code{nn.fit} is the object returned by \code{emfit}.
}
\references{
Rossell D. (2009) GaGa: a Parsimonious and Flexible Model for Differential Expression Analysis. Annals of Applied Statistics, 3, 1035-1051.

Yuan, M. and Kendziorski, C. (2006). A unified approach for simultaneous
gene clustering and differential expression identification.  Biometrics 62(4): 1089-1098.
}
\author{ David Rossell }
\seealso{ \code{\link{parest}} to estimate hyper-parameters and compute
  posterior probabilities after a GaGa or MiGaGa
  fit. \code{\link{findgenes}} to find differentially expressed
  genes. \code{\link{classpred}} to predict the group that a new sample
  belongs to. }
\examples{
library(gaga)
set.seed(10)
n <- 100; m <- c(6,6)
a0 <- 25.5; nu <- 0.109
balpha <- 1.183; nualpha <- 1683
probpat <- c(.95,.05)
xsim <- simGG(n,m,p.de=probpat[2],a0,nu,balpha,nualpha,equalcv=TRUE)
x <- exprs(xsim)

#Frequentist fit: EM algorithm to obtain MLE
groups <- pData(xsim)$group[c(-6,-12)]
patterns <- matrix(c(0,0,0,1),2,2)
colnames(patterns) <- c('group 1','group 2')
gg1 <- fitGG(x[,c(-6,-12)],groups,patterns=patterns,method='EM',trace=FALSE)  
gg1 <- parest(gg1,x=x[,c(-6,-12)],groups)
gg1

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
