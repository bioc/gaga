findgenes.gagafit <- function(gg.fit,x,groups,fdrmax=.05,parametric=TRUE,B=500) {
# Computes optimal terminal decision rules and expected terminal utility for several utility functions
# Input:
# - x: matrix or exprSet gene expression data
# - groups: vector indicating to which group does each column of x belong to
# - gg.fit: GaGa or MiGaGa fit, as returned by parest.gagafit
# - fdrmax: restriction on E(FDR). Ignored if util!='fnrstfdr'.
# - parametric: should the FDR be controlled parametrically or non-parametrically
# - B: number of permutations to estimate FDR non-parametrically (ignored if parametric==TRUE)
# - centers: ignored if parametric==TRUE. Indicates the number of clusters of z-scores (passed to kmeans). Note: centers==0 indicates to use no clusters (i.e. gene-specific bootstrap), while centers==1 indicates pooling residuals from all genes as in Storey's procedure.

# Output: list containing
# - efp: expected number of true positives
# - d: pattern to which each gene is assigned to
# - fdr: frequentist estimated FDR that is closest to fdrmax.
# - fdrpar: This is the target FDR used as input for Mueller's optimal Bayesian rule. If parametric==TRUE, this is equal to fdrmax. If parametric==FALSE, it's the Bayesian FDR needed to achieve frequentist estimated FDR=fdrmax.
# - fdrest: list with estimated frequentist FDR for each target Bayesian FDR
# - fnr: Bayesian FNR
# - power: Bayesian power as estimated by E(TP)/E(positives)
# - threshold: If util=='fnrstfdr', optimal threshold for posterior prob of pattern 0 (genes with prob < threshold are declared DE). Ignored otherwise.
# Note: Bayesian estimates are computed with model-based posterior probabilities. Frequentist estimates are obtained under repeated bootstrap sampling


centers <- 1; v <- gg.fit$pp
if (!is.matrix(v)) stop('gg.fit$pp must be a matrix containing posterior probabilities of each expression pattern')
cf <- as.double(2)

if (is(x, "exprSet") | is(x,"ExpressionSet")) {
  if (is.character(groups)) { groups <- as.factor(pData(data)[, groups]) }
  x <- exprs(x)
} else if (!is(x,"data.frame") & !is(x,"matrix")) { stop("x must be an exprSet, data.frame or matrix") }

groups <- as.integer(as.integer(as.factor(groups))-1); K <- as.integer(max(groups)+1)
if (length(groups) != ncol(x)) stop('groups must have length equal to the number of columns in x')
if (K==1) stop('At least two different groups must be specified')
patterns <- gg.fit$patterns
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
  znp <- .C("expected_fp",efp=fdrest,fdrseq,as.integer(length(fdrseq)),as.integer(B),as.integer(niter),as.double(t(zscore)),as.double(m),as.double(s),index,znclust,zclustsize,as.integer(nrow(x)),as.integer(ncol(x)),as.double(t(x)),as.integer(groups),as.integer(ncol(patterns)),as.double(a0),as.double(nu),as.double(balpha),as.double(nualpha),as.integer(gg.fit$equalcv),as.integer(length(probclus)),cluslist,as.double(t(probclus)),as.double(t(probpat)),as.integer(nrow(patterns)),as.integer(t(patterns)),ngrouppat,sumx,prodx,nobsx,sumxpred,prodxpred,nobsxpred,as.integer(gapprox))
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

return(list(efp=z$u,d=z$d,fdr=fdr,fdrpar=fdrpar,fdrest=fdrest,fnr=z$fnr,power=z$power,threshold=z$threshold))

}
