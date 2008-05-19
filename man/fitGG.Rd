\name{fitGG}
\alias{fitGG}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Fit GaGa hierarchical model }
\description{
 Fits GaGa or MiGaGa hierarchical models, either via a fully Bayesian
 approach or via maximum likelihood.
}
\usage{
fitGG(x, groups, patterns, equalcv = TRUE, nclust = 1, method = "SA", B, priorpar, parini, trace = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{\code{ExpressionSet}, \code{exprSet}, data frame or matrix
    containing the gene expression measurements used to fit the model.}
  \item{groups}{If \code{x} is of type \code{ExpressionSet} or
    \code{exprSet}, \code{groups} should be the name of the column
    in \code{pData(x)} with the groups that one wishes to compare. If
    \code{x} is a matrix or a data frame, \code{groups} should be a
    vector indicating to which group each column in x
    corresponds to.}
  \item{patterns}{Matrix indicating which groups are put together under
    each pattern, i.e. the hypotheses to consider for each
    gene. Defaults to two hypotheses: null hypothesis of all groups
    being equal and full alternative of all groups being different.}
  \item{equalcv}{\code{equalcv==TRUE} fits model assuming constant CV across groups. \code{equalcv==FALSE} compares cv as well as mean expression levels between groups}
  \item{nclust}{Number of clusters in the MiGaGa model. \code{nclust}
    corresponds to the GaGa model. }
  \item{method}{ \code{method=='MH'} fits a fully Bayesian model via
    Metropolis-Hastings posterior sampling. \code{method=='Gibbs'} does
    the same using Gibbs sampling. \code{method=='SA'} uses Simulated
    Annealing to find the posterior mode. \code{method=='EM'} finds
    maximum-likelihood estimates via the expectation-maximization
    algorithm, but this is currently only implemented for \code{nclust>1}. }
  \item{B}{ For \code{method=='MH'} and \code{method=='Gibbs'}, \code{B}
    is the number of MCMC iterations (defaults to 1000). For
    \code{method=='SA'}, \code{B} is the number of iterations in the
    Simulated Annealing scheme (defaults to 200). For
    \code{method=='EM'}, \code{B} is ignored. }
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
}
\details{
  An approximation is used to sample faster from the
  posterior distribution of the gamma shape parameters and to compute
  the normalization constants (needed to evaluate the likelihood). These
  approximations are implemented in \code{rcgamma} and \code{mcgamma}.

  The cooling scheme in \code{method=='SA'} uses a temperature equal to
  \code{1/log(1+i)}, where \code{i} is the iteration number.
}
\value{ An object of class \code{gagafit}, with components
  \item{parest }{Hyper-parameter estimates. Only returned if \code{method=='EBayes'}, for \code{method=='Bayes'} one must call the function \code{parest} after \code{fitGG}}
  \item{mcmc }{Object of class \code{mcmc} with posterior draws for hyper-parameters. Only returned if \code{method=='Bayes'}.}
  \item{lhood}{For \code{method=='Bayes'} it is the log-likelihood evaluated at each MCMC iteration. For \code{method=='EBayes'} it is the log-likelihood evaluated at the maximum.}
  \item{nclust}{Same as input argument.}
  \item{patterns}{Same as input argument, converted to object of class \code{gagahyp}.}
}
\references{ Rossell D. GaGa: a simple and 
flexible hierarchical model for microarray
data analysis. \url{http://rosselldavid.googlepages.com}. }
\author{ David Rossell }
\seealso{ \code{\link{parest}} to estimate hyper-parameters and compute
  posterior probabilities after a GaGa or MiGaGa
  fit. \code{\link{findgenes}} to find differentially expressed
  genes. \code{\link{classpred}} to predict the group that a new sample
  belongs to. }
\examples{
library(gaga)
library(EBarrays); data(gould)
x <- log(gould[1:100,-1])  #select 100 genes so that it runs very fast. Exclude 1st array
groups <- c(1,1,2,2,2,2,2,3,3) #group each array corresponds to
patterns <- rbind(rep(0,3),c(0,0,1),c(0,1,1),0:2) #4 hypothesis

#Frequentist fit: EM algorithm to obtain MLE
gg1 <- fitGG(x,groups,patterns,method='EM')
gg1 <- parest(gg1,x,groups)
gg1

#Bayesian fit: Metropolis-Hastings
gg.bay <- fitGG(x,groups,patterns,method='Gibbs',B=1000)
plot(gg.bay$mcmc)
gg.bay <- parest(gg.bay,x,groups,burnin=100)
gg.bay
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
