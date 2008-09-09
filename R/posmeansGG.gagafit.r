posmeansGG.gagafit <- function(gg.fit,x,groups,sel,underpattern) {
# Posterior means for each gene, assuming that pattern==underpattern holds for genes

gapprox <- TRUE
if (is(x, "exprSet") | is(x,"ExpressionSet")) {
  if (is.character(groups) && length(groups)==1) { groups <- as.factor(pData(x)[, groups]) }
  x <- exprs(x)
} else if (!is(x,"data.frame") & !is(x,"matrix")) { stop("x must be an exprSet, data.frame or matrix") } 

patterns <- gg.fit$patterns
v <- gg.fit$pp

groupsr <- groups2int(groups,patterns); K <- as.integer(max(groupsr)+1)

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
nclust <- as.integer(gg.fit$nclust); rho <- as.double(par$probclus)
sumx <- double(nrow(x)*sum(ngrouppat)); nobsx <- double(sum(ngrouppat))
if (gg.fit$equalcv) {
  prodx <- double(nrow(x))
} else {
  prodx <- double(nrow(x)*sum(ngrouppat))
}
gapprox <- as.integer(gapprox)

if (missing(underpattern)) {
  underpattern <- nrow(gg.fit$patterns)-1
} else if (underpattern>nrow(gg.fit$patterns)) {
  stop('The specified pattern number is not valid')
}
underpattern <- as.integer(underpattern)
cat(' Computing posterior means under expression pattern',underpattern,'...\n')

if (missing(sel)) sel <- 1:nrow(x)
if (is.logical(sel)) sel <- (1:nrow(x))[sel]
sel <- as.integer(sel-1) #in C indices start at 0
nsel <- length(sel)
posmeans <- double(nsel*K)

z <- .C("compute_sumxC",sumx=t(sumx),prodx=t(prodx),nobsx=nobsx,as.integer(gg.fit$equalcv),nsel,sel,as.integer(sum(ngrouppat)),ncol(x),as.double(t(x)),groupsr,K,as.integer(nrow(patterns)),as.integer(t(patterns)),ngrouppat,as.integer(1))
sumx <- matrix(z$sumx,nrow=nrow(x),byrow=TRUE); prodx <- matrix(z$prodx,nrow=nrow(x),byrow=TRUE)
nobsx <- z$nobsx

z <- .C("posmeans_ggC",posmeans=posmeans,underpattern,K,nsel,sel,alpha0,nu,balpha,nualpha,as.integer(gg.fit$equalcv),nclust,rho,as.double(t(v)),as.integer(nrow(patterns)),as.integer(t(patterns)),ngrouppat,t(sumx),t(prodx),nobsx,gapprox)

x <- data.frame(matrix(z$posmeans,nrow=nsel,byrow=TRUE))
for (i in 1:K) { names(x)[i] <- as.character(groups[groupsr==i-1][1]) }
return(x)

}
