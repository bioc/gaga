\name{seqBoundariesGrid}
\alias{seqBoundariesGrid}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Evaluate expected utility for parametric sequential stopping boundaries. }
\description{
  Estimate the expected utility for sequential boundaries
  parameterized by (b0,b1). Expected utility is estimated on a grid of
  (b0,b1) values based on a forward simulation output such as that
  generated by the function \code{forwsimDiffExpr}.
}
\usage{
seqBoundariesGrid(b0, b1, forwsim, samplingCost, powmin = 0, f = "linear", ineq = "less")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{b0}{ Vector with b0 values. Expected utility is evaluated for a
    grid defined by all combinations of (b0,b1) values. }
  \item{b1}{ Vector with b1 values. }
  \item{forwsim}{ \code{data.frame} with forward simulation output, such
  as that returned by the function \code{forwsimDiffExpr}. It must have
  columns named \code{simid}, \code{time}, \code{u}, \code{fdr},
  \code{fnr}, \code{power} and \code{summary}. See
  \code{forwsimDiffExpr} for details on the meaning of each column.}
  \item{samplingCost}{ Cost of obtaining one more data batch, in terms
    of the number of new truly differentially expressed discoveries that
  would make it worthwhile to obtain one more data batch.}
  \item{powmin}{ Constraint on power. Optimization chooses the optimal
    \code{b0}, \code{b1} satisfying power>=powermin (if such \code{b0},\code{b1} exists). }
  \item{f}{ Parametric form for the stopping boundary. Currently only
    'linear' and 'invsqrt' are implemented. For 'linear', the boundary
    is \code{b0+b1*time}. For 'invsqrt', the boundary is
    \code{b0+b1/sqrt(time)}, where time is the sample size measured as
    number of batches. }
  \item{ineq}{ For \code{ineq=='less'} the trial stops when \code{summary} is below
  the stopping boundary. This is appropriate whenever \code{summary}
  measures the potential benefit of obtaining one more data batch. For
  \code{ineq=='greater'} the trial stops when \code{summary} is above
  the stopping boundary. This is approapriate whenever \code{summary}
  measures the potential costs of obtaining one more data batch.}
}
\details{
  Intuitively, the goal is to stop collecting new data when the expected
  benefit of obtaining one more data batch is small, i.e. below a
  certain boundary. We consider two simple parametric forms for such a
  boundary (linear and inverse square root), which allows to easily evaluate
  the expected utility for each boundary within a grid of parameter
  values.
  The optimal boundary is defined by the parameter values achieving the
  largest expected utility, restricted to parameter values with an
  estimated power greater or equal than \code{powmin}.
  Here power is defined as the expected number of true discoveries
  divided by the expected number of differentially expressed entities.

  The routine evaluates the expected utility, as well as expected FDR,
  FNR, power and sample size for each specified boundary, and also reports
  the optimal boundary.
}
\value{
  A list with two components:
  \item{opt }{Vector with optimal stopping boundary (\code{b}), estimated expected
    utility (\code{u}), false discovery rate (\code{fdr}), false
    negative rate (\code{fnr}), power (\code{power}) and the expected
    sample size measured as the number of batches (\code{time}).}
  \item{grid }{\code{data.frame} with all evaluated boundaries (columns \code{b0}
    and \code{b1}) and their respective estimated expected utility,
    false discovery rate, false negative rate, power and expected sample
  size (measured as the number of batches).}
}
\references{
  Rossell D., Mueller P. Sequential sample sizes for
  high-throughput hypothesis testing
  experiments. \url{http://sites.google.com/site/rosselldavid/home}.
  
  Rossell D. GaGa: a simple and flexible hierarchical model for microarray
   data analysis. Annals of Applied Statistics, 2009, 3, 1035-1051.
}
\author{ David Rossell. }
\seealso{ \code{\link{forwsimDiffExpr}} }
\keyword{ htest }
\keyword{ design }% __ONLY ONE__ keyword per line
