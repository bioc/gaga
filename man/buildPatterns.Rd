\name{buildPatterns}
\alias{buildPatterns}
\title{
Build a matrix with all possible patterns given a number of groups where
samples may belong to.
}
\description{
Creates a matrix indicating which groups are put together under each
pattern. The number of possible patterns increases very fast as the
number of groups increases. This function provides an easy way to
compute all possible patterns.
The output of this function is usually used for the \code{patterns}
parameter of the \code{lmFit} function.
}
\usage{buildPatterns(groups)}
\arguments{
\item{groups}{Character containing the names of the groups at which
  samples may belong to.
  If the output of the function is going to be used in \code{fitGG} it
  must match the group levels specified in the \code{groups} parameter
  that will be given to \code{fitGG}.}
}
\examples{
buildPatterns(groups=c('GroupControl','GroupA','GroupB'))
}
\keyword{logic}
