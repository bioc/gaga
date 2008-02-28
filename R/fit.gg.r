fit.gg <- function(x,groups,patterns,nclust=1,method='Bayes',B=1000,priorpar,parini,trace=TRUE) {
# Fits GaGa or MiGaGa model via MCMC or empirical Bayes
# Input:
# - x: matrix with gene expression measurements for all groups
# - groups: vector of length ncol(x) indicating what group does each column in x correspond to
# - patterns: matrix indicating which groups are put together under each pattern, i.e. the hypotheses to consider for each gene. Defaults to two hypotheses: null hypothesis of all groups being equal and full alternative of all groups being different
# - nclust: number of clusters in hierarchical prior
# - method: 'Bayes' fits model via Gibbs' MCMC posterior simulation, 'EBayes' fits model via empirical Bayes
# - B: number of MCMC samples to obtain when method=='Bayes'. Ignored for method=='EBayes'.
# - priorpar: list with prior parameter values (a.alpha0,b.alpha0,a.nu,b.nu,a.balpha,b.balpha,a.nualpha,b.nualpha,p.probclus,p.probpat).
# - parini: list with initial values for hyper-parameters a0, nu, balpha, nualpha, probclus, probpat
# Output: an object of class gagafit, with components
# - parest: parameter estimates. Only returned if method=='EBayes', for method=='Bayes' one must call the function parest after fit.gg
# - mcmc: posterior draws for hyper-parameters. Only returned if method=='Bayes'.
# - lhood: for method=='Bayes' it is the log-likelihood evaluated at each MCMC iteration. For method=='EBayes' it is the log-likelihood evaluated at the maximum
# - nclust: number of clusters
# - patterns: same as input argument

gapprox <- TRUE
if (is(x, "ExpressionSet")) {
  if (is.character(groups)) { groups <- as.factor(pData(x)[, groups]) }
  x <- exprs(x)
} else if (!is(x,"data.frame") & !is(x,"matrix")) { stop("x must be an ExpressionSet, data.frame or matrix") }
groups <- as.integer(as.integer(as.factor(groups))-1); K <- as.integer(max(groups)+1)
if (ncol(x)!=length(groups)) stop('length(groups) must be equal to the number of columns in x')
if (missing(patterns)) patterns <- rbind(rep(0,K),0:(K-1))
if (length(table(groups))!=ncol(patterns)) stop('patterns must have the same number of columns as the number of distinct groups')
if (ncol(patterns)!=K) stop('patterns must have number of columns equal to the number of distinct elements in groups')
for (i in 1:nrow(patterns)) { patterns[i,] <- as.integer(as.integer(as.factor(patterns[i,]))-1) }
class(patterns) <- 'gagahyp'
ngrouppat <- as.integer(apply(patterns,1,'max')+1)
nclust <- as.integer(nclust)
npat <- as.integer(nrow(patterns))

# Initialize parameters by method of moments
if (nclust==1) probclusini <- 1
if (missing(parini)) {
  if (trace) cat('Initializing parameters...')
  aest <- rowMeans(x)^2/apply(x,1,'var'); lest <- 1/rowMeans(x)
  sel <- (aest<quantile(aest,probs=.99)) & (lest<quantile(lest,probs=.99))
  aest <- aest[sel]; lest <- lest[sel]
  balphaini <- as.double(mean(aest)^2/var(aest)); nualphaini <- as.double(mean(aest))
  if (nclust==1) {
    a0ini <- as.double(mean(lest)^2/var(lest)); nuini <- as.double(mean(lest))
    probclusini <- as.double(1)
  } else {
    clusini <- kmeans(x=lest,centers=nclust)$cluster
    a0ini <- as.double(tapply(lest,clusini,'mean')^2/tapply(lest,clusini,'var'))
    nuini <- as.double(tapply(lest,clusini,'mean'))
    probclusini <- as.double(table(clusini)/length(clusini))
  }
  probpatini <- rep(1/K,K)
  if (trace) cat(' Done.\n')
} else {
  if (is.null(parini$a0) | is.null(parini$nu) | is.null(parini$balpha) | is.null(parini$nualpha) | is.null(parini$probpat)) stop('some components of parini are empty')
  a0ini <- parini$a0; nuini <- parini$nu; balphaini <- parini$balpha; nualphaini <- parini$nualpha
  if (nclust==1) { probclusini <- 1 } else { if (is.null(parini$probclus)) { stop('component probclus of parini is empty') } else { probclusini <- parini$probclus/sum(parini$probclus) } }
  probpatini <- parini$probpat/sum(parini$probpat)
}


if (method=='EBayes') {

  if (nclust>1) stop('nclust>1 is not currently implemented for empirical Bayes method')
  if (trace) cat('Starting EM algorithm...\n')
  z <- fitfreq.gg(x,groups,patterns,nclust=1,a0ini,nuini,balphaini,nualphaini,probpatini,iter.max=100,trace=trace)
  parest <- c(alpha0=z$alpha0,nu=z$nu,balpha=z$balpha,nualpha=z$nualpha,probclus=1,probpat.1=z$probpat[1],probpat.2=z$probpat[2])
  gg.fit <- list(parest=parest,mcmc=as.mcmc(NA),lhood=z$lhood,nclust=nclust,patterns=patterns,method=method)
  class(gg.fit) <- 'gagafit'
  return(gg.fit)

} else if (method=='Bayes') {

  if (missing(priorpar)) {
    a.alpha0 <-  .0016; b.alpha0 <-  .0001; a.nu <-  .016; b.nu <-  .0016
    a.balpha <-  .004; b.balpha <- .001; a.nualpha <- .004; b.nualpha <- 20
    p.probclus <- as.double(rep(.1,nclust))
    p.probpat <- as.double(rep(.1,npat))
  } else {
    if (is.null(priorpar$a.alpha0) | (is.null(priorpar$b.alpha0)) | (is.null(priorpar$a.nu)) | (is.null(priorpar$b.nu)) | (is.null(priorpar$a.balpha)) | (is.null(priorpar$b.balpha)) | (is.null(priorpar$a.nualpha)) | (is.null(priorpar$b.nualpha))) stop('Some components of priorpar have not been specified')
    a.alpha0 <- priorpar$a.alpha0; b.alpha0 <- priorpar$b.alpha0
    a.nu <- priorpar$a.nu; b.nu <- priorpar$b.nu
    a.balpha <- priorpar$a.balpha; b.balpha <- priorpar$b.balpha
    a.nualpha <- priorpar$a.nualpha; b.nualpha <- priorpar$b.nualpha
    if (nclust>1) { p.probclus <- 1 } else { if (is.null(priorpar$p.probclus)) stop('component p.probclus of priorpar has not been specified') else p.probclus <- priorpar$p.probclus }
    p.probpat <- priorpar$p.probpat
  }

# Call MCMC sampling routine
  if (trace) cat('MCMC sampling...\n')
  balpha <- nualpha <- double(B); alpha0 <- nu <- probclus <- double(B*nclust); prob <- double(B*npat)
  lhood <- double(B); trace <- as.integer(trace)

  z <- .C("fit_ggC",alpha0=alpha0,nu=nu,balpha=balpha,nualpha=nualpha,probclus=probclus,prob=prob,lhood=lhood,as.integer(B),as.double(a.alpha0),as.double(b.alpha0),as.double(a.nu),as.double(b.nu),as.double(a.balpha),as.double(b.balpha),as.double(a.nualpha),as.double(b.nualpha),as.double(p.probclus),as.double(p.probpat),a0ini,nuini,balphaini,nualphaini,probclusini,probpatini,as.integer(nrow(x)),as.integer(ncol(x)),as.double(t(x)),as.integer(groups),K,nclust,npat,as.integer(t(patterns)),ngrouppat,as.integer(gapprox),trace)
  if (trace) cat('Done.\n')

  gg.fit <- list(parest=NA,mcmc=as.mcmc(data.frame(alpha0=matrix(z$alpha0,nrow=B,ncol=nclust,byrow=TRUE),nu=matrix(z$nu,nrow=B,ncol=nclust,byrow=TRUE),balpha=z$balpha,nualpha=z$nualpha,probclus=matrix(z$probclus,nrow=B,ncol=nclust,byrow=TRUE),probpat=matrix(z$prob,nrow=B,byrow=TRUE))),lhood=z$lhood,nclust=nclust,patterns=patterns,method=method)
  class(gg.fit) <- 'gagafit'
  return(gg.fit)

}

}
