\name{print.gagaclus}
\alias{print.gagaclus}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Print an object of class gagaclus }
\description{
  Prints an object of class \code{gagaclus}, which contains the result
  of clustering genes into expression patterns.
}
\usage{
print.gagaclus(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ Object of type \code{gagaclus}. }
  \item{\dots}{Other arguments to be passed on to the generic print function.}
}
\details{
}
\value{
  Displays the expression patterns and the number of genes classified into
  each of them.
}
\references{ Rossell D. GaGa: a simple and 
flexible hierarchical model for microarray
data analysis. \url{http://rosselldavid.googlepages.com}. }
\author{ David Rossell }
\seealso{ \code{\link{fitGG}}, \code{\link{geneclus}} }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ print }
