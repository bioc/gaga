\name{print.gagahyp}
\alias{print.gagahyp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Print an object of class gagahyp }
\description{
  Prints an object of class \code{gagahyp}, which contains information
  on the hypotheses (expression patterns) from a GaGa or MiGaGa model.
}
\usage{
\method{print}{gagahyp}(x, probpat=NA, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ Object of type \code{gagahyp}. }
  \item{probpat}{ Vector with either estimated probabilities of each
    hypothesis, or with number of genes classified into each expression pattern. }
  \item{\dots}{Other arguments to be passed on to the generic print function.}
}
\value{
Prints hypotheses. When available, also displays estimated proportion of
genes following each expression pattern or the number of genes classified into each
expression pattern.
}
\references{ Rossell D. GaGa: a simple and 
flexible hierarchical model for microarray
data analysis. \url{http://rosselldavid.googlepages.com}. }
\author{ David Rossell }
\seealso{ \code{\link{fitGG}}, \code{\link{geneclus}} }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ print }
