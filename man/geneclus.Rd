\name{geneclus}
\alias{geneclus}
\alias{geneclus.gagafit}
\title{ Cluster genes into expression patterns. }
\description{
  Performs supervised gene clustering. Clusters genes into the expression pattern with highest posterior probability,
  according to a GaGa or MiGaGa fit.
}
\usage{
geneclus(gg.fit, method='posprob')
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gg.fit}{GaGa or MiGaGa fit (object of type \code{gagafit}, as
    returned by \code{fitGG}). }
  \item{method}{For \code{method==1} samples are assigned to pattern
    with highest posterior probability, and for \code{method==1} to the
    pattern with highest 
    likelihood (e.g. assuming equal a priori prob for all patterns)}
}
\details{
  Each gene is assigned to the pattern with highest posterior
  probability.  
  This is similar to routine \code{findgenes}, which also assigns genes to
  the pattern with highest posterior probability, although
  \code{findgenes} applies an FDR-based correction i.e. tends to assign
  more genes to the null pattern of no differential expression.
}
\value{
List with components:
  \item{d }{Vector indicating the pattern that each gene is assigned
to.}
  \item{posprob }{Vector with posterior probabilities of the assigned patterns.}
}
\references{ Rossell D. GaGa: a simple and 
flexible hierarchical model for microarray
data analysis. \url{http://rosselldavid.googlepages.com}. }
\author{ David Rossell }
\seealso{ \code{\link{fitGG}}, \code{\link{parest}} }
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
#dclus <- geneclus(ggfit)  #not use FDR correction
#dfdr <- findgenes(ggfit,xsim$x[,c(-6,-12)],groups,fdrmax=.05,parametric=TRUE) #use FDR correction
#table(dfdr$d,dclus$d) #compare results
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }
\keyword{ models }
