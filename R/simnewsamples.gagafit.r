simnewsamples.gagafit <- function(gg.fit,groupsnew,sel,x,groups) {
# Simulates parameter values and new observations for a select subset of genes and the group indicated by groups from a Gamma/Gamma model
# - groupsnew: vector indicating the groups that each new observation will belong to
# - sel: numeric vector indicating the indexes of the genes we want to draw new samples for (defaults to all genes). If a logical vector is indicated, it is converted to (1:nrow(x))[sel]
# - x: matrix with observations used to fit the model
# - groups: vector of length ncol(x) indicating what group does each column in x correspond to
# - gg.fit: GaGa fit, as returned by parest.gagafit
# Output
# - xnew: matrix of nsel rows and nsamples cols with obs drawn from the predictive for groups given by groupsnew.
# - dnew: matrix of nsel rows and nsamples cols with expression pattern indicators drawn from the posterior.
# - anew: matrix of nsel rows and nsamples cols with alpha parameters from the posterior for groups given by groupsnew.
# - lnew: matrix of nsel rows and nsamples cols with lambda parameters from the posterior for groups given by groupsnew.

gapprox <- TRUE
if (is(x,"ExpressionSet")) {
  if (is.character(groups)) { groups <- as.factor(pData(data)[, groups]) }
  x <- exprs(x)
} else if (!is(x,"data.frame") & !is(x,"matrix")) { stop("x must be an ExpressionSet, data.frame or matrix") }

patterns <- gg.fit$patterns
v <- gg.fit$pp
if ((max(groupsnew)>max(groups)) | (min(groupsnew)<min(groups))) stop('Groups indicated in groupsnew do not match with those indicated in groups')
groupstemp <- factor(groups)
groupsnew <- factor(groupsnew,levels=levels(groupstemp))
groups <- as.integer(as.integer(groupstemp)-1); K <- as.integer(max(groups)+1)
groupsnew <- as.integer(as.integer(groupsnew)-1)
if (ncol(x)!=length(groups)) stop('length(groups) must be equal to the number of columns in x')
if (!is.matrix(v)) stop('Argument v must be a matrix')
if (ncol(v)!=nrow(patterns)) stop('Argument v must have as many columns as rows has patterns')
if (nrow(x)!=nrow(v)) stop('Arguments x and v must have the same number of rows')

if (ncol(patterns)!=K) stop('patterns must have number of columns equal to the number of distinct elements in groups')
for (i in 1:nrow(patterns)) { patterns[i,] <- as.integer(as.integer(as.factor(patterns[i,]))-1) }
ngrouppat <- as.integer(apply(patterns,1,'max')+1)
par <- getpar(gg.fit)
alpha0 <- as.double(par$a0); nu <- as.double(par$nu); balpha <- as.double(par$balpha)
nualpha <- as.double(par$nualpha)
nclust <- as.integer(length(gg.fit$a0)); rho <- as.double(par$probclus)
sumx <- prodx <- double(nrow(x)*sum(ngrouppat)); nobsx <- double(sum(ngrouppat))
gapprox <- as.integer(gapprox)

if (missing(sel)) sel <- 1:nrow(x)
if (is.logical(sel)) sel <- (1:nrow(x))[sel]
sel <- as.integer(sel-1) #in C indices start at 0
nsel <- length(sel); nsamples <- length(groupsnew)
xnew <- anew <- lnew <- double(nsel*nsamples); dnew <- integer(nsel*nsamples)

z <- .C("compute_sumxC",sumx=t(sumx),prodx=t(prodx),nobsx=nobsx,nsel,sel,as.integer(sum(ngrouppat)),ncol(x),as.double(t(x)),groups,K,as.integer(nrow(patterns)),as.integer(t(patterns)),ngrouppat,as.integer(1))
sumx <- matrix(z$sumx,nrow=nrow(x),byrow=TRUE); prodx <- matrix(z$prodx,nrow=nrow(x),byrow=TRUE)
nobsx <- z$nobsx

z <- .C("simnewsamples_ggC",xnew=xnew,dnew=dnew,anew=anew,lnew=lnew,nsamples,groupsnew,K,nsel,sel,alpha0,nu,balpha,nualpha,nclust,rho,as.double(t(v)),as.integer(nrow(patterns)),as.integer(t(patterns)),ngrouppat,t(sumx),t(prodx),nobsx,gapprox)

return(list(xnew=matrix(z$xnew,nrow=nsel,byrow=TRUE),dnew=matrix(z$dnew,nrow=nsel,byrow=TRUE),anew=matrix(z$anew,nrow=nsel,byrow=TRUE),lnew=matrix(z$lnew,nrow=nsel,byrow=TRUE)))

}
