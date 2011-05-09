\name{print.gagafit}
\alias{print.gagafit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Print an object of class gagafit }
\description{
  Prints an object of class \code{gagafit}, as returned by \code{fitGG}
  or \code{parest}. Provides general information and hyper-parameter
  estimates, if available.
}
\usage{
\method{print}{gagafit}(x,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ Object of type \code{gagafit}, as returned by \code{fitGG}
    or \code{parest}. }
  \item{\dots}{Other arguments to be passed on to the generic print function.}
}
\details{
  \code{fitGG} does not create a complete \code{gagafit} object. The
  complete object is returned by \code{parest}, which computes the
  posterior probabilities of differential expression and obtain
  hyper-parameter estimates (these are only provided by \code{fitGG}
  when the option \code{method='EBayes'} is used).
}
\value{
Prints number of genes, hypotheses, details about the model fitting and
hyper-parameter estimates (when available).
}
\references{ Rossell D. GaGa: a simple and 
flexible hierarchical model for microarray
data analysis. \url{http://rosselldavid.googlepages.com}. }
\author{ David Rossell }
\seealso{ \code{\link{fitGG}}, \code{\link{parest}} }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ print }
