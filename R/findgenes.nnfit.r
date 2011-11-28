findgenes.nnfit <- function(fit,x,groups,fdrmax=.05,parametric=TRUE,B=500) {
#Check input
v <- fit$pp
if (!is.matrix(v)) stop('fit$pp must be a matrix containing posterior probabilities of each expression pattern')
if (parametric==FALSE) stop("parametric==FALSE not implemented for Normal-Normal model")
cf <- as.double(2)
if (is(x, "exprSet") | is(x,"ExpressionSet")) {
  if (is.character(groups) && length(groups)==1) { groups <- as.factor(pData(x)[, groups]) }
  x <- exprs(x)
} else if (!is(x,"data.frame") & !is(x,"matrix")) { stop("x must be an exprSet, data.frame or matrix") }
patterns <- fit$patterns
groupsr <- groups2int(groups,patterns); K <- as.integer(max(groupsr)+1)
if (length(groups) != ncol(x)) stop('groups must have length equal to the number of columns in x')
if (K==1) stop('At least two different groups must be specified')
if (ncol(patterns)!=K) stop('patterns must have number of columns equal to the number of distinct elements in groups')
#Find genes
nsel <- nrow(v); sel <- as.integer((1:nsel)-1)
util <- as.integer(1)
u <- double(1); d <- integer(nrow(v)); fdr <- fnr <- power <- double(1); threshold <- double(1)
z <- .C("utgene_parC",u=u,d=d,fdr=fdr,fnr=fnr,power=power,threshold=threshold,util=util,as.double(cf),nsel,sel,as.double(t(v)),as.integer(ncol(v)),as.double(fdrmax))
fdr <- fdrest <- fdrpar <- z$fdr
return(list(truePos=z$u,d=z$d,fdr=fdr,fdrpar=fdrpar,fdrest=fdrest,fnr=z$fnr,power=z$power,threshold=z$threshold))
}
