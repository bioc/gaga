\name{sim.gg}
\alias{sim.gg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Prior predictive simulation }
\description{
  Simulates parameters and data from the prior-predictive of GaGa or
  MiGaGa model with 2 groups, fixing the hyper-parameters.
}
\usage{
sim.gg(n, m, p.de=.1, a0, nu, balpha, nualpha, probclus=1, a=NA, l=NA, useal=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{ Number of genes. }
  \item{m}{ Vector indicating number of observations to be simulated for
  each group.}
  \item{p.de}{ Probability that a gene is differentially expressed. }
  \item{a0, nu}{ Mean expression for each gene is generated from
  \code{1/rgamma(a0,a0/nu)} if \code{probclus} is of length 1, and from a
  mixture if \code{length(probclus)>1}.}
  \item{balpha, nualpha}{ Shape parameter for each gene is generated from \code{rgamma(balpha,balpha/nualpha)}.}
  \item{probclus}{Vector with the probability of each component in the
  mixture. Set to 1 for the GaGa model.}
  \item{a, l}{ Optionally, if \code{useal==TRUE} the parameter values are
  not generated, only the data is generated. \code{a} is a matrix with the shape parameters
  of each gene and group and \code{l} is a matrix with the mean expressions. }
  \item{useal}{ For \code{useal==TRUE} the parameter values specified in
  \code{a} and \code{l} are used, instead of being generated.}
}
\value{
  A list with the following components:
  \item{x}{Matrix with simulated observations (group 1 in columns
    \code{1:m[1]}, group 2 in columns \code{(m[1]+1):(m[1]+m[2])} etc.)}
  \item{a}{Matrix with shape parameters used to generate x.}
  \item{l}{Matrix with mean expression used to generate x.}
}
\details{
  The shape parameters are actually drawn from a gamma approximation to
  their posterior distribution. The function \code{rcgamma} implements
  this approximation.
}
\references{ Rossell D. GaGa: a simple and 
flexible hierarchical model for microarray
data analysis. \url{http://rosselldavid.googlepages.com}. }
\author{ David Rossell }
\note{ Currently, the routine only implements prior predictive
  simulation for the 2 hypothesis case.
}
\seealso{ \code{\link{simnewsamples}} to simulate from the posterior
  predictive, \code{\link{checkfit}} for graphical posterior predictive checks. }
\examples{
#Not run. Example from the help manual
#library(gaga)
#set.seed(10)
#n <- 100; m <- c(6,6)
#a0 <- 25.5; nu <- 0.109
#balpha <- 1.183; nualpha <- 1683
#probpat <- c(.95,.05)
#xsim <- sim.gg(n,m,p.de=probpat[2],a0,nu,balpha,nualpha)
#
#plot(density(xsim$x),main='')
#plot(xsim$l,xsim$a,ylab='Shape',xlab='Mean')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ distribution }
\keyword{ models }
