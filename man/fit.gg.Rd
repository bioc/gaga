\name{fit.gg}
\alias{fit.gg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Fit GaGa hierarchical model}
\description{
 Fits GaGa or MiGaGa hierarchical models, either via a fully Bayesian
 approach or via maximum likelihood.
}
\usage{
fit.gg(x, groups, patterns, nclust = 1, method = "Bayes", B = 1000, priorpar, parini, trace = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{\code{ExpressionSet}, data frame or matrix
    containing the gene expression measurements used to fit the model.}
  \item{groups}{If \code{x} is of type \code{ExpressionSet},
    \code{groups} should be the name of the column in \code{pData(x)}
    with the groups that one wishes to compare. If \code{x} is a matrix
    or a data frame, \code{groups} should be a vector indicating to which
    group each column in x corresponds to.}
  \item{patterns}{Matrix indicating which groups are put together under each pattern, i.e. the hypotheses to consider for each gene. Defaults to two hypotheses: null hypothesis of all groups being equal and full alternative of all groups being different.}
  \item{nclust}{Number of clusters in the MiGaGa model. \code{nclust}
    corresponds to the GaGa model.}
  \item{method}{\code{method=='Bayes'} fits a fully Bayesian model via
    MCMC posterior sampling. \code{method=='EBayes'} finds
    maximum-likelihood estimates via the expectation-maximization
    algorithm. For \code{nclust>1} only \code{Bayes} is currently implemented.}
  \item{B}{Number of MCMC iterations. Ignored if \code{method=='EBayes'}.}
  \item{priorpar}{List with prior parameter values. It must have
    components \code{a.alpha0,b.alpha0,a.nu,b.nu,a.balpha,b.balpha,a.nualpha,b.nualpha,p.probclus}
    and \code{p.probpat}. If missing they are set to non-informative
    values that are usually reasonable for RMA and GCRMA normalized data.}
  \item{parini}{List with components \code{a0}, \code{nu},
    \code{balpha}, \code{nualpha}, \code{probclus} and \code{probpat}
    indicating the starting values for the hyper-parameters. If not
    specified, a method of moments estimate is used.}
  \item{trace}{For \code{trace==TRUE} the progress of the model fitting
  routine is printed.}
}
\details{
  The Bayesian fit uses an approximation to sample faster from the
  posterior distribution of the gamma shape parameters. This
  approximation is implemented in \code{rcgamma}.
}
\value{An object of class \code{gagafit}, with components
  \item{parest }{Hyper-parameter estimates. Only returned if \code{method=='EBayes'}, for \code{method=='Bayes'} one must call the function \code{parest} after \code{fit.gg}}
  \item{mcmc }{Object of class \code{mcmc} with posterior draws for hyper-parameters. Only returned if \code{method=='Bayes'}.}
  \item{lhood}{For \code{method=='Bayes'} it is the log-likelihood evaluated at each MCMC iteration. For \code{method=='EBayes'} it is the log-likelihood evaluated at the maximum.}
  \item{nclust}{Same as input argument.}
  \item{patterns}{Same as input argument, converted to object of class \code{gagahyp}.}
}
\references{Rossell D. GaGa: a simple and 
flexible hierarchical model for microarray
data analysis. \url{http://rosselldavid.googlepages.com}.}
\author{David Rossell}
\seealso{\code{\link{parest}} to estimate hyper-parameters and compute
  posterior probabilities after a GaGa or MiGaGa
  fit. \code{\link{findgenes}} to find differentially expressed
  genes. \code{\link{classpred}} to predict the group that a new sample
  belongs to.}
\examples{
#Not run
#library(EBarrays); data(gould)
#x <- log(exprs(gould)[,-1])  #exclude 1st array
#groups <- pData(gould)[-1,1]
#patterns <- rbind(rep(0,3),c(0,0,1),c(0,1,1),0:2) #4 hypothesis
#gg <- fit.gg(x,groups,patterns,method='EBayes')
#gg <- parest(gg,x,groups)
#gg
#
#gg.bay <- fit.gg(x,groups,patterns,method='Bayes',B=1000)
#plot(gg.bay$mcmc)
#gg.bay <- parest(gg.bay,x,groups,burnin=100)
#gg.bay
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
