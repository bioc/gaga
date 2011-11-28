\name{findgenes}
\alias{findgenes}
\alias{findgenes.gagafit}
\alias{findgenes.nnfit}
\title{ Find differentially expressed genes after GaGa or Normal-Normal fit. }
\description{
  Obtains a list of differentially expressed genes using the posterior
  probabilities from a GaGa, MiGaGa or Normal-Normal fit. For \code{parametric==TRUE} the procedure
  controls the Bayesian FDR below \code{fdrmax}. For
  \code{parametric==FALSE} it controls the estimated frequentist FDR
  (only available for GaGa).
}
\usage{
findgenes(fit, x, groups, fdrmax=.05, parametric=TRUE, B=500)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{Either GaGa/MiGaGa fit (object of class \code{gagafit}, as
    returned by \code{fitGG}) or Normal-Normal fit (class \code{nnfit},
    as returned by \code{fitNN}). }
  \item{x}{\code{ExpressionSet}, \code{exprSet}, data frame or matrix
    containing the gene expression measurements used to fit the model.}
  \item{groups}{If \code{x} is of type \code{ExpressionSet} or
    \code{exprSet}, \code{groups} should be the name of the column
    in \code{pData(x)} with the groups that one wishes to compare. If
    \code{x} is a matrix or a data frame, \code{groups} should be a
    vector indicating to which group each column in x
    corresponds to.}
  \item{fdrmax}{Upper bound on FDR.}.
  \item{parametric}{Set to \code{TRUE} to control the posterior expected
    FDR below \code{fdrmax}. Set to
    \code{FALSE} to estimate the frequentist FDR non-parametrically
    (only available when fit is of class \code{gagafit}).}
  \item{B}{Number of boostrap samples to estimate FDR non-parametrically (ignored if \code{parametric==TRUE})}
}
\details{
  The Bayes rule to minimize posterior expected FNR subject to posterior
  expected FDR
  \code{<=fdrmax} declares differentially expressed all genes with
  posterior probability of being equally expressed below a certain
  threshold. The value of the threshold is computed exactly for
  \code{parametric==TRUE}, FDR being defined in a Bayesian sense. For
  \code{parametric==FALSE} the FDR is defined in a frequentist sense.
}
\value{
List with components:
  \item{truePos }{Expected number of true positives.}
  \item{d }{Vector indicating the pattern that each gene is assigned
to.}
  \item{fdr }{Frequentist estimated FDR that is closest to fdrmax.}
  \item{fdrpar}{Bayesian FDR. If \code{parametric==TRUE}, this is equal to \code{fdrmax}. If \code{parametric==FALSE}, it's the Bayesian FDR needed to achieve frequentist estimated FDR=\code{fdrmax}.}
  \item{fdrest}{Data frame with estimated frequentist FDR for each target Bayesian FDR}
  \item{fnr}{Bayesian FNR}
  \item{power}{Bayesian power as estimated by expected number of true
    positives divided by the expected number of differentially expressed
  genes}
  \item{threshold}{Optimal threshold for posterior probability of equal expression (genes with probability < \code{threshold} are declared DE)}
}
\references{
Rossell D. (2009) GaGa: a Parsimonious and Flexible Model for Differential Expression Analysis. Annals of Applied Statistics, 3, 1035-1051.

Yuan, M. and Kendziorski, C. (2006). A unified approach for simultaneous
gene clustering and differential expression identification.  Biometrics 62(4): 1089-1098.

Muller P, Parmigiani G, Robert C, Rousseau J. (2004) Journal of the American
Statistical Association, 99(468): 990-1001.
}
\author{ David Rossell }
\seealso{ \code{\link{fitGG}}, \code{\link{fitNN}}, \code{\link{parest}} }
\examples{
#Not run. Example from the help manual
#library(gaga)
#set.seed(10)
#n <- 100; m <- c(6,6)
#a0 <- 25.5; nu <- 0.109
#balpha <- 1.183; nualpha <- 1683
#probpat <- c(.95,.05)
#xsim <- simGG(n,m,p.de=probpat[2],a0,nu,balpha,nualpha)
#
#ggfit <- fitGG(xsim$x[,c(-6,-12)],groups,patterns=patterns,nclust=1)
#ggfit <- parest(ggfit,x=xsim$x[,c(-6,-12)],groups,burnin=100,alpha=.05)
#
#d <- findgenes(ggfit,xsim$x[,c(-6,-12)],groups,fdrmax=.05,parametric=TRUE)
#dtrue <- (xsim$l[,1]!=xsim$l[,2])
#table(d$d,dtrue)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }
\keyword{ models }
