findgenes.nnfit <- function(fit,x,groups,fdrmax=.05,parametric=TRUE,B=500) {
#Check input
if (!is.matrix(fit$pp)) stop('fit$pp must be a matrix containing posterior probabilities of each expression pattern')
if (parametric==FALSE) stop("parametric==FALSE not implemented for Normal-Normal model")
cf <- as.double(2)
#Find genes
nsel <- nrow(fit$pp); sel <- as.integer((1:nsel)-1)
util <- as.integer(1)
u <- double(1); d <- integer(nrow(fit$pp)); fdr <- fnr <- power <- double(1); threshold <- double(1)
z <- .C("utgene_parC",u=u,d=d,fdr=fdr,fnr=fnr,power=power,threshold=threshold,util=util,as.double(cf),nsel,sel,as.double(t(fit$pp)),as.integer(ncol(fit$pp)),as.double(fdrmax))
fdr <- fdrest <- fdrpar <- z$fdr
return(list(truePos=z$u,d=z$d,fdr=fdr,fdrpar=fdrpar,fdrest=fdrest,fnr=z$fnr,power=z$power,threshold=z$threshold))
}
