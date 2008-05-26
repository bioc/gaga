ppGG <- function(x,groups,a0,nu,balpha,nualpha,equalcv,probclus,probpat,patterns) {
# Posterior probabilities of DE in GaGa and MiGaGa models given data and hyper-param estimate

gapprox <- TRUE
if (is(x, "exprSet") | is(x,"ExpressionSet")) {
  if (is.character(groups)) { groups <- as.factor(pData(data)[, groups]) }
  x <- exprs(x)
} else if (!is(x,"data.frame") & !is(x,"matrix")) { stop("x must be an exprSet, data.frame or matrix") }

groupsr <- groups2int(groups,patterns)
K <- as.integer(max(groupsr)+1)
if (ncol(x)!=length(groups)) stop('length(groups) must be equal to the number of columns in x')
if (missing(a0)) stop('a0 must be specified')
if (missing(nu)) stop('nu must be specified')
if (missing(balpha)) stop('balpha must be specified')
if (missing(nualpha)) stop('nualpha must be specified')
if (!is.vector(probclus)) stop('probclus must be a vector')
if (!is.vector(probpat)) stop('probpat must be a vector')
if ((length(a0)!=length(nu)) || (length(a0)!=length(probclus))) stop('a0,nu and probclus must have the same length')
if (length(balpha)>1 || length(nualpha)>1) stop('balpha and nualpha must be vectors of length 1')

if (ncol(patterns)!=K) stop('patterns must have number of columns equal to the number of distinct elements in groups')
for (i in 1:nrow(patterns)) { patterns[i,] <- as.integer(as.integer(as.factor(patterns[i,]))-1) }
ngrouppat <- as.integer(apply(patterns,1,'max')+1)
v <- double(nrow(x)*nrow(patterns)); lhood <- double(1)
usesumx <- as.integer(0)
sumx <- double(nrow(x)*sum(ngrouppat)); nobsx <- double(sum(ngrouppat))
sumxpred <- double(nrow(x)*sum(ngrouppat));  nobsxpred <- double(sum(ngrouppat));
if (as.logical(equalcv)) {
  prodx <- double(nrow(x));
  prodxpred <- double(nrow(x));
} else {
  prodx <- double(nrow(x)*sum(ngrouppat));
  prodxpred <- double(nrow(x)*sum(ngrouppat));
}
nsel <- nrow(x); sel <- as.integer(0:(nsel-1))
cluslist <- as.integer(c((0:(length(probclus)-1)),-1))
z <- .C("pp_ggC",v=v,lhood=lhood,nsel,sel,as.integer(ncol(x)),as.double(t(x)),groupsr,as.integer(ncol(patterns)),as.double(a0),as.double(nu),as.double(balpha),as.double(nualpha),as.integer(equalcv),as.integer(length(probclus)),cluslist,as.double(t(probclus)),as.double(t(probpat)),as.integer(nrow(patterns)),as.integer(t(patterns)),ngrouppat,sumx,prodx,nobsx,sumxpred,prodxpred,nobsxpred,usesumx,as.integer(gapprox))
v <- matrix(z$v,nrow=nrow(x),byrow=TRUE)

return(list(pp=v,lhood=z$lhood))

}
