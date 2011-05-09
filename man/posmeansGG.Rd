\name{posmeansGG}
\alias{posmeansGG}
\alias{posmeansGG.gagafit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Gene-specific posterior means }
\description{
  Computes posterior means for the gene expression levels using a GaGa or MiGaGa model.
}
\usage{
posmeansGG(gg.fit, x, groups, sel, underpattern)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gg.fit}{GaGa or MiGaGa fit (object of type \code{gagafit}, as returned by \code{fitGG}). }
  \item{x}{\code{ExpressionSet}, \code{exprSet}, data frame or matrix
    containing the gene expression measurements used to fit the model.}
  \item{groups}{If \code{x} is of type \code{ExpressionSet} or
    \code{exprSet}, \code{groups} should be the name of the column
    in \code{pData(x)} with the groups that one wishes to compare. If
    \code{x} is a matrix or a data frame, \code{groups} should be a
    vector indicating to which group each column in x
    corresponds to.}
  \item{sel}{Numeric vector with the indexes of the genes we want to
    draw new samples for (defaults to all genes). If a logical vector is
    indicated, it is converted to \code{(1:nrow(x))[sel]}.}
  \item{underpattern}{Expression pattern assumed to be true (defaults to
    last pattern in \code{gg.fit$patterns}). Posterior
    means are computed under this pattern. For example, if only the null
    pattern that all groups are equal and the full alternative that all
    groups are different are considered, \code{underpattern=1} returns
    the posterior means under the assumption that groups are different
    from each other (\code{underpattern=0} returns the same mean for all
  groups).}
}
\details{
  The posterior distribution of the mean parameters actually depends on
  the gene-specific shape parameter(s), which is unknown. To speed up
  computations, a gamma approximation to the shape parameter posterior
  is used (see \code{rcgamma} for details) and the shape parameter is
  fixed to its mode a posteriori.
}
\value{
  Matrix with mean expression values a posteriori, for each selected
  gene and each group. Genes are in rows and groups in columns.
}
\references{ Rossell D. GaGa: a simple and 
flexible hierarchical model for microarray
data analysis. \url{http://rosselldavid.googlepages.com}. }
\author{ David Rossell }
\seealso{ \code{\link{fitGG}} for fitting GaGa and MiGaGa models,
  \code{\link{parest}} for computing posterior probabilities of
  each expression pattern. }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ distribution }
\keyword{ models }
