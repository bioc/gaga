\name{plotForwSim}
\alias{plotForwSim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plot forward simulation trajectories }
\description{
  Produces plot to visualize simulated trajectories of the summary
  statistic as a function of the number of additional data batches,
  i.e. the output produced by \code{forwsimDiffExpr}.
}
\usage{
plotForwSim(fs,xlab="Number of additional batches per group",ylab="Expected increase in True Positives",col='gray',lty=1,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fs}{Forward simulation results, as output by \code{forwsimDiffExpr}.}
  \item{xlab}{x-axis label}
  \item{ylab}{y-axis label}
  \item{col}{line color}
  \item{lty}{line type}
  \item{\dots}{ Other arguments to be passed to \code{plot} }
}
\value{
 Produces a plot.
}
\references{   Rossell D., Mueller P. Sequential sample sizes for
  high-throughput hypothesis testing
  experiments. \url{http://sites.google.com/site/rosselldavid/home}.
 }
\author{ David Rossell }
\seealso{ \code{\link{forwsimDiffExpr}} }
\keyword{ plots }
