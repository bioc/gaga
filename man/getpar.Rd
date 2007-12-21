\name{getpar}
\alias{getpar}
\alias{getpar.gagafit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Extract hyper-parameter estimates from a gagafit object }
\description{
  Extracts the hyper-parameter estimates from a \code{gagafit} object
  and puts them in a list.
}
\usage{
getpar(gg.fit)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gg.fit}{ Object of class \code{gagafit}, as returned by \code{parest}. }
}
\details{
  This routine simply evaluates the component \code{gg.fit$parest} from
  a \code{gagafit} object, which causes an error if this component
  is not available. This routine is used internally by a number of other
  routines.
}
\value{
  A list with components:
  \item{a0 }{Estimated value of hyper-parameter \code{a0}}
  \item{nu }{Estimated value of hyper-parameter \code{nu}}
  \item{balpha }{Estimated value of hyper-parameter \code{balpha}}
  \item{nualpha }{Estimated value of hyper-parameter \code{nualpha}}
  \item{probclus }{Estimated cluster probabilities}
  \item{probpat }{Estimated prior probability of each expression pattern}
}
\references{ Rossell D. GaGa: a simple and 
flexible hierarchical model for microarray
data analysis. \url{http://rosselldavid.googlepages.com}. }
\author{ David Rossell }
\seealso{ \code{\link{fit.gg}}, \code{\link{parest}} }
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
