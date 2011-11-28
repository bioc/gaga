\name{getpar}
\alias{getpar}
\alias{getpar.gagafit}
\alias{getpar.nnfit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Extract hyper-parameter estimates from a gagafit or nnfit object }
\description{
  Extracts the hyper-parameter estimates from a \code{gagafit} or
  \code{nnfit} object
  and puts them in a list.
}
\usage{
getpar(fit)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{ Object of class \code{gagafit} or \code{nnfit}. }
}
\details{
  This routine simply evaluates the component \code{parest} from
  a \code{gagafit} or \code{nnfit} object, which causes an error if this component
  is not available. This routine is used internally by a number of other
  routines.
}
\value{
  For \code{gagafit} objects, a list with components:
  \item{a0 }{Estimated value of hyper-parameter \code{a0}}
  \item{nu }{Estimated value of hyper-parameter \code{nu}}
  \item{balpha }{Estimated value of hyper-parameter \code{balpha}}
  \item{nualpha }{Estimated value of hyper-parameter \code{nualpha}}
  \item{probclus }{Estimated cluster probabilities}
  \item{probpat }{Estimated prior probability of each expression
    pattern}

  For \code{nnfit} objects, a vector with elements \code{mu0},
  \code{tau02}, \code{v0}, \code{sigma02}, \code{probclus} and
  \code{probpat}. These are the hierarchical N(mu0,tau0^2) *
  IG(.5*v0,.5*v0*sigma0^2) prior parameter estimates.
}
\references{
  Rossell D. (2009) GaGa: a Parsimonious and Flexible Model for Differential Expression Analysis. Annals of Applied Statistics, 3, 1035-1051.

Yuan, M. and Kendziorski, C. (2006). A unified approach for simultaneous
gene clustering and differential expression identification.  Biometrics 62(4): 1089-1098.
   }
\author{ David Rossell }
\seealso{ \code{\link{fitGG}}, \code{\link{fitNN}}, \code{\link{parest}} }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
