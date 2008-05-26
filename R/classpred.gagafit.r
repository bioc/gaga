classpred.gagafit <- function(gg.fit,xnew,x,groups,prgroups,ngene=100) {
# Classify a new sample (array) into the group with highest posterior probability

gapprox <- TRUE
sel <- (1:nrow(x))[order(gg.fit$pp[,1])][1:ngene]
if (!is.vector(xnew)) stop('xnew must be a vector')
if (!is.numeric(sel)) stop('sel must contain numerical indexes')
if (is(x, "exprSet") | is(x,"ExpressionSet")) {
  if (is.character(groups)) { groups <- as.factor(pData(data)[, groups]) }
  x <- exprs(x)
} else if (!is(x,"data.frame") & !is(x,"matrix")) { stop("x must be an exprSet, data.frame or matrix") }
if (ncol(x)!=length(groups)) stop('Argument groups must have length equal to number of columns in argument x')
par <- getpar(gg.fit)
a0 <- as.double(par$a0); nu <- as.double(par$nu); balpha <- as.double(par$balpha)
nualpha <- as.double(par$nualpha)
probclus <- as.double(par$probclus); probpat <- as.double(par$probpat)
equalcv <- as.integer(gg.fit$equalcv)
nclust <- as.integer(length(probclus))
patterns <- gg.fit$patterns
if (nrow(patterns)!=length(probpat)) stop('Argument patterns must be equal to the length of gg.fit@probEst')

groupsr <- groups2int(groups,patterns); K <- as.integer(max(groupsr)+1)
if (missing(prgroups)) prgroups <- rep(1/K,K)
prgroups <- as.double(prgroups)
if (ncol(patterns)!=K) stop('patterns must have number of columns equal to the number of distinct elements in groups')
npat <- as.integer(nrow(patterns))
for (i in 1:nrow(patterns)) { patterns[i,] <- as.integer(as.integer(as.factor(patterns[i,]))-1) }
ngrouppat <- as.integer(apply(patterns,1,'max')+1)
sumx <- double(nrow(x)*sum(ngrouppat)); nobsx <- double(sum(ngrouppat))
if (gg.fit$equalcv) {
  prodx <- double(nrow(x))
} else {
  prodx <- double(nrow(x)*sum(ngrouppat))
}
usesumx <- as.integer(0); gapprox <- as.integer(gapprox)

xnew <- as.double(xnew[sel])
sel <- sel-1  #in C vectors start at 0
d <- integer(1); posgroups <- double(K)

z <- .C("sampleclas_ggC",d=d,posgroups=posgroups,xnew,as.integer(length(sel)),as.integer(sel),as.integer(nrow(x)),as.integer(ncol(x)),as.double(t(x)),groupsr,K,prgroups,probclus,probpat,a0,nu,balpha,nualpha,equalcv,nclust,npat,as.integer(t(patterns)),ngrouppat,sumx,prodx,nobsx,usesumx,gapprox)

return(list(d=z$d+1,posgroups=z$posgroups))

}
