powclasspred.gagafit <- function(gg.fit,x,groups,prgroups,v0thre=1,ngene=100,B=100) {
# Estimate expected probability that a future sample is correctly classified.

gapprox <- TRUE
patterns <- gg.fit$patterns
if (is(x, "exprSet") | is(x,"ExpressionSet")) {
  if (is.character(groups) && length(groups)==1) { groups <- as.factor(pData(x)[, groups]) }
  x <- exprs(x)
} else if (!is(x,"data.frame") & !is(x,"matrix")) { stop("x must be an exprSet, data.frame or matrix") }
if (ncol(x)!=length(groups)) stop('Argument groups must have length equal to number of columns in argument x')
par <- getpar(gg.fit)
a0 <- as.double(par$a0); nu <- as.double(par$nu); balpha <- as.double(par$balpha)
nualpha <- as.double(par$nualpha); probclus <- as.double(par$probclus); probpat <- as.double(par$probpat); nclust <- as.integer(gg.fit$nclust)
if (nrow(patterns)!=length(probpat)) stop('Argument patterns must be equal to the length of gg.fit$prob')
if (missing(genelimit)) { genelimit <- nrow(x); }

genelimit <- as.integer(genelimit); v0thre <- as.double(v0thre)
usesel <- as.integer(0); nsel <- as.integer(0); sel <- integer(nrow(x))

npat <- as.integer(nrow(patterns))
groupsr <- groups2int(groups,patterns); K <- as.integer(max(groupsr)+1)
if (missing(prgroups)) prgroups <- rep(1/K,K)
if (ncol(patterns)!=K) stop('patterns must have number of columns equal to the number of distinct elements in groups')
for (i in 1:nrow(patterns)) { patterns[i,] <- as.integer(as.integer(as.factor(patterns[i,]))-1) }
ngrouppat <- as.integer(apply(patterns,1,'max')+1)
ncolsumx <- as.integer(sum(ngrouppat))
sumx <- double(nrow(x)*ncolsumx); nobsx <- double(ncolsumx)
if (gg.fit$equalcv) {
  prodx <- double(nrow(x))
} else {
  prodx <- double(nrow(x)*sum(ngrouppat))
}
usesumx <- as.integer(0); gapprox <- as.integer(gapprox)

ccall <- seccall <- double(1); ccgroup <- double(K); ngroup <- integer(K); preceps <- as.double(0)

v <- ppGG(x=x,groups=groups,a0=a0,nu=nu,balpha=balpha,nualpha=nualpha,equalcv=gg.fit$equalcv,probclus=probclus,probpat=probpat,patterns=patterns)$pp

z <- .C("utsample_ggC",ccall=ccall,seccall=seccall,ccgroup=ccgroup,ngroup=ngroup,as.integer(B),preceps,genelimit,v0thre,nsel,sel,usesel,as.integer(nrow(x)),as.integer(ncol(x)),as.double(t(x)),as.integer(groupsr),as.double(t(v)),K,as.double(prgroups),a0,nu,balpha,nualpha,as.integer(gg.fit$equalcv),nclust,probclus,probpat,npat,as.integer(t(patterns)),ngrouppat,ncolsumx,sumx,prodx,nobsx,usesumx,gapprox)

ccgroup <- z$ccgroup/z$ngroup
return(list(ccall=z$ccall,seccall=z$seccall,ccgroup=ccgroup,segroup=ccgroup*(1-ccgroup)/sqrt(z$ngroup)))

}
