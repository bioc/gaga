\name{dcgamma}
\alias{dcgamma}
\alias{rcgamma}
\alias{mcgamma}
\title{ Approximate gamma shape distribution }
\description{
  \code{dcgamma} approximates density of a gamma shape distribution with
  a gamma density. \code{rcgamma} obtains random draws from the
  approximation. \code{mcgamma} computes approximated mean, variance and
  normalization constant.
}
\usage{
dcgamma(x, a, b, c, d, r, s, newton = TRUE)
rcgamma(n, a, b, c, d, r, s, newton = TRUE)
mcgamma(a, b, c, d, r, s, newton = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ Vector indicating the values at which to evaluate the
    density. }
  \item{n}{ Number of random draws to obtain.}
  \item{a,b,c,d,r,s}{ Parameter values. }
  \item{newton}{ Set to \code{TRUE} to try to locate the mode by taking
    a few Newton-Raphson steps. }  
}
\details{
  The density of a gamma shape distribution is given by
  \code{C(a,b,c,d,r,s) (gamma(a*x+d)/gamma(x)^a)
  (x/(r+s*x))^{a*x+d} x^{b-d-1} exp(-x*c)}
  for \code{x>=0}, and 0 otherwise, where \code{C()} is the normalization constant.
  The gamma approximation is
  \code{Ga(a/2+b-1/2,c+a*log(s/a))}. The approximate normalization constant is
  obtained by taking the ratio of the exact density and the
  approximation at the maximum, as described in Rossell (2007).
}
\value{
\code{dcgamma} returns a vector with approximate density.
\code{rcgamma} returns a vector with draws from the approximating gamma.
\code{mcgamma} returns a list with components:
  \item{m }{Approximate mean}
  \item{v }{Approximate variance}
  \item{normk }{Approximate normalization constant}
}
\references{ Rossell D. GaGa: a simple and 
flexible hierarchical model for microarray
data analysis. \url{http://rosselldavid.googlepages.com}. }
\author{ David Rossell }
\note{ For general values of the parameters the gamma approximation may
  be poor. In such a case one could use this function to obtain draws
  from the proposal distribution in a Metropolis-Hastings step. 
}
\seealso{ \code{\link{dgamma}}, \code{\link{rgamma}} }
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ distribution }
