\name{simnewsamples}
\alias{simnewsamples}
\alias{simnewsamples.gagafit}
\alias{simnewsamples.nnfit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Posterior predictive simulation }
\description{
  Posterior and posterior predictive simulation for GaGa/MiGaGa and
  Normal-Normal models.
}
\usage{
simnewsamples(fit, groupsnew, sel, x, groups)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{Either GaGa or MiGaGa fit (object of type \code{gagafit},
    as returned by \code{fitGG}) or Normal-Normal fit (type \code{nnfit}
    returned by \code{fitNN}). }
  \item{groupsnew}{ Vector indicating the group that each new sample
    should belong to. \code{length(groupsnew)} is the number of new
    samples that will be generated. }
  \item{sel}{Numeric vector with the indexes of the genes we want to
    draw new samples for (defaults to all genes). If a logical vector is
    indicated, it is converted to \code{(1:nrow(x))[sel]}. For the
    Normal-Normal model this argument is ignored.}
  \item{x}{\code{ExpressionSet}, \code{exprSet}, data frame or matrix
    containing the gene expression measurements used to fit the model.}
  \item{groups}{If \code{x} is of type \code{ExpressionSet} or
    \code{exprSet}, \code{groups} should be the name of the column
    in \code{pData(x)} with the groups that one wishes to compare. If
    \code{x} is a matrix or a data frame, \code{groups} should be a
    vector indicating to which group each column in x
    corresponds to.}
}
\details{
  For GaGa/MiGaGa models, the shape parameters are actually drawn from a gamma approximation to
  their posterior distribution. The function \code{rcgamma} implements
  this approximation.

  In order to be consistent with the LNNGV model implemented in emfit (package
  EBarrays), for the Normal-Normal model the variance is drawn from an inverse
  gamma approximation to its marginal posterior (obtained by plugging in
  the group means, see EBarrays vignette for details).
}
\value{
  Object of class 'ExpressionSet'. Expression values can be accessed via
  \code{exprs(object)} and the parameter values used to generate the
  expression values can be accessed via \code{fData(object)}.
}
\references{
Rossell D. (2009) GaGa: a Parsimonious and Flexible Model for Differential Expression Analysis. Annals of Applied Statistics, 3, 1035-1051.

Yuan, M. and Kendziorski, C. (2006). A unified approach for simultaneous
gene clustering and differential expression identification.  Biometrics 62(4): 1089-1098.
}
\author{ David Rossell }
\seealso{ \code{\link{checkfit}} for posterior predictive plot,
  \code{\link{simGG}} for prior predictive simulation. }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ distribution }
\keyword{ models }
