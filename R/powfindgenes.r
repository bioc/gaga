# Estimates expected predictive terminal utility after drawing batchSize more observations per group from the predictive distrib
powfindgenes <- function(fit,x,groups,batchSize=1,fdrmax=.05,genelimit,v0thre=1,B=1000,mc.cores=1) { UseMethod("powfindgenes") }

powfindgenes.nnfit <- function(fit,x,groups,batchSize=1,fdrmax=.05,genelimit,v0thre=1,B=1000,mc.cores=1) {
#Format input
if (is(x, "exprSet") | is(x,"ExpressionSet")) {
  if (is.character(groups) && length(groups)==1) { groups <- as.factor(pData(x)[, groups]) }
  x <- exprs(x)
} else if (!is(x,"data.frame") & !is(x,"matrix")) { stop("x must be an exprSet, data.frame or matrix") }
groupsnew <- rep(colnames(fit$patterns),each=batchSize)
probpat <- getpar(fit); probpat <- probpat[grep('probpat',names(probpat))]
#Function to simulate + post prob + post expected true positives
f <- function(simid) {
  xnew <- simnewsamples(fit,groupsnew=groupsnew,x=x,groups=groups)
  fitall <- updateNNfit(fit,x=cbind(x,exprs(xnew)),groups=c(groups,groupsnew))
  class(fitall) <- 'nnfit'
  findgenes(fitall,fdrmax=fdrmax)$truePos
}
#Do simulations
if (mc.cores>1) {
  if ('multicore' %in% loadedNamespaces()) {
    ans <- multicore::mclapply(1:B, f, mc.cores=mc.cores, mc.preschedule=TRUE)
  } else stop('multicore library has not been loaded!')
} else {
  ans <- lapply(1:B, f)
}
ans <- unlist(ans)
return(list(m=mean(ans),s=sd(ans)/sqrt(B)))
}

powfindgenes.gagafit <- function(fit,x,groups,batchSize=1,fdrmax=.05,genelimit,v0thre=1,B=1000,mc.cores=1) {
gapprox <- TRUE
patterns <- fit$patterns
if (is(x, "exprSet") | is(x,"ExpressionSet")) {
  if (is.character(groups) && length(groups)==1) { groups <- as.factor(pData(x)[, groups]) }
  x <- exprs(x)
} else if (!is(x,"data.frame") & !is(x,"matrix")) { stop("x must be an exprSet, data.frame or matrix") }
if (ncol(x)!=length(groups)) stop('Argument groups must have length equal to number of columns in argument x')
v <- fit$pp
if (!is.matrix(v)) stop('Argument v must be a matrix')
if (nrow(x)!=nrow(v)) stop('Arguments x and v must have the same number of rows')
if (nrow(patterns)!=ncol(v)) stop('Argument patterns must have number of rows equal to the number of columns in v')
if (missing(genelimit)) { genelimit <- nrow(x); }

genelimit <- as.integer(genelimit); v0thre <- as.double(v0thre)
usesel <- as.integer(0); nsel <- as.integer(0); sel <- integer(nrow(x))

par <- getpar(fit)
alpha0 <- as.double(par$a0); nu <- as.double(par$nu); balpha <- as.double(par$balpha)
nualpha <- as.double(par$nualpha)
probclus <- as.double(par$probclus); probpat <- as.double(par$probpat)
if (ncol(v)!=length(probpat)) stop('Argument prob must have length equal to number of columns in v')

m <- s <- double(1)
util <- as.integer(1); cf <- double(2)
npat <- as.integer(length(probpat))
groupsr <- groups2int(groups,patterns); K <- as.integer(max(groupsr)+1)
if (ncol(patterns)!=K) stop('patterns must have number of columns equal to the number of distinct elements in groups')
for (i in 1:nrow(patterns)) { patterns[i,] <- as.integer(as.integer(as.factor(patterns[i,]))-1) }
ngrouppat <- as.integer(apply(patterns,1,'max')+1)

sumx <- double(nrow(x)*sum(ngrouppat)); nobsx <- double(sum(ngrouppat))
if (fit$equalcv) {
  prodx <- double(nrow(x))
} else {
  prodx <- double(nrow(x)*sum(ngrouppat))
}
preceps <- as.double(0); usesumx <- as.integer(0); gapprox <- as.integer(gapprox)

z <- .C("utgene_predC",m=m,s=s,as.integer(batchSize),as.integer(B),preceps,util,cf,genelimit,v0thre,nsel,sel,usesel,as.integer(nrow(x)),as.integer(ncol(x)),as.double(t(x)),as.integer(groupsr),K,as.double(t(v)),alpha0,nu,balpha,nualpha,as.integer(fit$equalcv),as.integer(length(probclus)),probclus,probpat,npat,as.integer(t(patterns)),ngrouppat,as.double(fdrmax),sumx,prodx,nobsx,usesumx,gapprox)

return(list(m=z$m,s=z$s))

}
