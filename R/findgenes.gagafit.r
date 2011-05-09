findgenes.gagafit <- function(gg.fit,x,groups,fdrmax=.05,parametric=TRUE,B=500) {

centers <- 1; v <- gg.fit$pp
if (!is.matrix(v)) stop('gg.fit$pp must be a matrix containing posterior probabilities of each expression pattern')
cf <- as.double(2)

if (is(x, "exprSet") | is(x,"ExpressionSet")) {
  if (is.character(groups) && length(groups)==1) { groups <- as.factor(pData(x)[, groups]) }
  x <- exprs(x)
} else if (!is(x,"data.frame") & !is(x,"matrix")) { stop("x must be an exprSet, data.frame or matrix") }

patterns <- gg.fit$patterns
groupsr <- groups2int(groups,patterns); K <- as.integer(max(groupsr)+1)
if (length(groups) != ncol(x)) stop('groups must have length equal to the number of columns in x')
if (K==1) stop('At least two different groups must be specified')
if (ncol(patterns)!=K) stop('patterns must have number of columns equal to the number of distinct elements in groups')
for (i in 1:nrow(patterns)) { patterns[i,] <- as.integer(as.integer(as.factor(patterns[i,]))-1) }
ngrouppat <- as.integer(apply(patterns,1,'max')+1)
par <- getpar(gg.fit)
a0 <- as.double(par$a0); nu <- as.double(par$nu); balpha <- as.double(par$balpha); nualpha <- as.double(par$nualpha); probclus <- as.double(par$probclus); probpat <- as.double(par$probpat)
cluslist <- as.integer(c((0:(length(probclus)-1)),-1))
if (B<10) { warning('B was set to less than 10, too small a number of permutations. Increased to B=10'); B <- 10 }

sumx <- double(nrow(x)*sum(ngrouppat)); nobsx <- double(sum(ngrouppat))
sumxpred <- double(nrow(x)*sum(ngrouppat)); nobsxpred <- double(sum(ngrouppat))
if (gg.fit$equalcv) {
  prodx <- prodxpred <- double(nrow(x))
} else {
  prodx <- prodxpred <- double(nrow(x)*sum(ngrouppat))
}
gapprox <- 1

nsel <- nrow(v); sel <- as.integer((1:nsel)-1)
util <- as.integer(1)
u <- double(1); d <- integer(nrow(v)); fdr <- fnr <- power <- double(1); threshold <- double(1)
z <- .C("utgene_parC",u=u,d=d,fdr=fdr,fnr=fnr,power=power,threshold=threshold,util=util,as.double(cf),nsel,sel,as.double(t(v)),as.integer(ncol(v)),as.double(fdrmax))
fdr <- fdrest <- fdrpar <- z$fdr

if (parametric==FALSE) {
  fdrseq <- as.double(seq(fdrmax/1000,min(fdrmax*2,1),length=1000))
  fdrest <- double(length(fdrseq))
  cat("Finding clusters of z-scores for bootstrap... ")
  m <- rowMeans(x); s <- sqrt((rowMeans(x^2)-rowMeans(x)^2)*ncol(x)/(ncol(x)-1))
  zscore <- (x-m)/s
  if (centers>1) {
    nquant <- min(10,ncol(x))
    qx <- t(matrix(unlist(apply(zscore,1,'quantile',probs=seq(0,1,length=nquant))),nrow=nquant))
    zcluster <- kmeans(qx,centers=centers,iter.max=50)
    zclustsize <- as.integer(table(zcluster$cluster))
    index <- as.integer(order(zcluster$cluster)-1)
  } else if ((centers==1) | (centers==0)) {
    zclustsize <- nrow(x)
    index <- as.integer(0:(nrow(x)-1))
  } else { stop('centers must be an integer >=0') }
  znclust <- as.integer(centers); niter <- 10

  cat("Done\nStarting",B,"bootstrap iterations...\n")
  znp <- .C("expected_fp",efp=fdrest,fdrseq,as.integer(length(fdrseq)),as.integer(B),as.integer(niter),as.double(t(zscore)),as.double(m),as.double(s),index,znclust,zclustsize,as.integer(nrow(x)),as.integer(ncol(x)),as.double(t(x)),as.integer(groupsr),as.integer(ncol(patterns)),as.double(a0),as.double(nu),as.double(balpha),as.double(nualpha),as.integer(gg.fit$equalcv),as.integer(length(probclus)),cluslist,as.double(t(probclus)),as.double(t(probpat)),as.integer(nrow(patterns)),as.integer(t(patterns)),ngrouppat,sumx,prodx,nobsx,sumxpred,prodxpred,nobsxpred,as.integer(gapprox))
  fdrest <- znp$efp
  fdrpar <- fdrseq[abs(fdrest-fdrmax)==min(abs(fdrest-fdrmax))][1]
  fdr <- fdrest[abs(fdrest-fdrmax)==min(abs(fdrest-fdrmax))][1]
  if (fdr>1.1*fdrmax) {
    warning('estimated FDR too large. Try decreasing fdrmax')
    #fdrpar <- 0
  }
  z <- .C("utgene_parC",u=u,d=d,fdr=fdr,fnr=fnr,power=power,threshold=threshold,util=util,as.double(cf),nsel,sel,as.double(t(v)),as.integer(ncol(v)),as.double(fdrpar))
  fdrest <- data.frame(fdrseq=fdrseq,fdrest=fdrest)
}

return(list(truePos=z$u,d=z$d,fdr=fdr,fdrpar=fdrpar,fdrest=fdrest,fnr=z$fnr,power=z$power,threshold=z$threshold))

}
