\name{parest}
\alias{parest}
\alias{parest.gagafit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Parameter estimates and posterior probabilities of differential expression for GaGa and MiGaGa model }
\description{
  Obtains parameter estimates and posterior probabilities of
  differential expression after a GaGa or MiGaGa model has been fit with
  the function \code{fit.gg}.
}
\usage{
parest(gg.fit, x, groups, burnin, alpha=.05)
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
  \item{burnin}{Number of MCMC samples to discard. Ignored if
    \code{gg.fit} was fit with the option \code{method=='EBayes'}.}
  \item{alpha}{If \code{gg.fit} was fit with the option
    \code{method=='Bayes'}, \code{parest} also computes \code{1-alpha}
    posterior credibility intervals.}
}
\details{
  If \code{gg.fit} was fit via MCMC posterior sampling (option
  \code{method=='Bayes'}), \code{parest} discards the first
  \code{burnin} iterations and uses the rest to obtain point estimates
  and credibility intervals for the hyper-parameters.
  To compute posterior probabilities of differential expression the hyper-parameters are fixed to
  their estimated value, i.e. not averaged over MCMC iterations.
}
\value{
  An object of class \code{gagafit}, with components:
  \item{parest }{Hyper-parameter estimates.}
  \item{mcmc }{Object of class \code{mcmc} with posterior draws for hyper-parameters. Only returned if \code{method=='Bayes'}.}
  \item{lhood}{For \code{method=='Bayes'} it is the posterior mean of the log-likelihood. For \code{method=='EBayes'} it is the log-likelihood evaluated at the maximum.}
  \item{nclust}{Number of clusters.}
  \item{patterns}{Object of class \code{gagahyp} indicating which
    hypotheses (expression patterns) were tested.}
  \item{pp}{Matrix with posterior probabilities of differential
    expression for each gene. Genes are in rows and expression patterns
    are in columns (e.g. for 2 hypotheses, 1st column is the probability of
    the null hypothesis and 2nd column for the alternative).}
}
\references{ Rossell D. GaGa: a simple and 
flexible hierarchical model for microarray
data analysis. \url{http://rosselldavid.googlepages.com}. }
\author{ David Rossell }
\seealso{ \code{\link{fit.gg}} to fit a GaGa or MiGaGa model,
  \code{\link{findgenes}} to find differentially expressed genes and
  \code{\link{classpred}} to do class prediction. }
\examples{
#Not run
#library(EBarrays); data(gould)
#x <- log(exprs(gould)[,-1])  #exclude 1st array
#groups <- pData(gould)[-1,1]
#patterns <- rbind(rep(0,3),c(0,0,1),c(0,1,1),0:2) #4 hypothesis
#gg <- fit.gg(x,groups,patterns,method='EBayes')
#gg
#gg <- parest(gg,x,groups)
#gg
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
