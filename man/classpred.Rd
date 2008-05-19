\name{classpred}
\alias{classpred}
\alias{classpred.gagafit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Predict the class that a new sample belongs to. }
\description{
  Computes the posterior probability that a new sample belongs to each
  group and classifies it into the group with highest probability.
}
\usage{
classpred(gg.fit, xnew, x, groups, prgroups, ngene=100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gg.fit}{GaGa or MiGaGa fit (object of type \code{gagafit}, as returned by \code{fitGG}). }
  \item{xnew}{Expression levels of the sample to be classified. Only the subset of the genes indicated by ngene is used.}
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
  \item{ngene}{Number of genes to use to build the classifier. Genes
    with smaller probability of being equally expressed are selected first.}
}
\details{
  The classifier weights each gene according to the posterior
  probability that it is differentially expressed. Hence, adding genes
  that are unlikely to be differentially expressed does not affect the
  performance of the classifier, but it does increase the computational
  cost.
  All computations are performed by fixing the hyper-parameters to their
  estimated value (posterior mean if model was fit with
  \code{method=='Bayes'} or maximum likelihood estimate is model was fit
  with \code{method=='EBayes'}).
}
\value{
  List with the following elements:
  \item{d }{Numeric value indicating the group that the new sample is
    classified into, i.e. where the maximum in \code{posgroups} is.} 
  \item{posgroups }{Vector giving the posterior probability that the
    \code{xnew} belongs to each of the groups.}
}
\references{ Rossell D. GaGa: a simple and 
flexible hierarchical model for microarray
data analysis. \url{http://rosselldavid.googlepages.com}. }
\author{ David Rossell }
\seealso{ \code{\link{fitGG}}, \code{\link{parest}} }
\examples{
#Not run. Example from the help manual
#library(gaga)
#set.seed(10)
#n <- 100; m <- c(6,6)
#a0 <- 25.5; nu <- 0.109
#balpha <- 1.183; nualpha <- 1683
#probpat <- c(.95,.05)
#xsim <- simGG(n,m,p.de=probpat[2],a0,nu,balpha,nualpha)
#
#ggfit <- fitGG(xsim$x[,c(-6,-12)],groups,patterns=patterns,nclust=1)
#ggfit <- parest(ggfit,x=xsim$x[,c(-6,-12)],groups,burnin=100,alpha=.05)
#
#pred1 <- classpred(ggfit,xnew=xsim$x[,6],x=xsim$x[,c(-6,-12)],groups)
#pred2 <- classpred(ggfit,xnew=xsim$x[,12],x=xsim$x[,c(-6,-12)],groups)
#pred1
#pred2
}
\keyword{ htest }
\keyword{ models }
