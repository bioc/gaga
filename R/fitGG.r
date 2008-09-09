fitGG <- function(x,groups,patterns,equalcv=TRUE,nclust=1,method='quickEM',B,priorpar,parini,trace=TRUE) {

#Input processing: check errors, format and set missing parameters to default
if (missing(B)) { if (method=='SA') { B <- 200 } else if (method=='MH' | method=='Gibbs') { B <- 1000 } else { B <- 10 } }
gapprox <- TRUE
if (is(x, "exprSet") | is(x, "ExpressionSet")) {
  if (is.character(groups) && length(groups)==1) { groups <- as.factor(pData(x)[, groups]) }
  x <- exprs(x)
} else if (!is(x,"data.frame") & !is(x,"matrix")) { stop("x must be an ExpressionSet, exprSet, data.frame or matrix") }
if (min(x)<0) stop("x can only have positive values")
if (sum(is.na(x))>0) stop("x cannot have any NA values")
if (ncol(x)!=length(groups)) stop('length(groups) must be equal to the number of columns in x')
if ((method=='quickEM') && (nrow(x)>10^4)) x <- x[sample(1:nrow(x),10^4),]  #limit nb genes to 10,000 for speed

K <- length(unique.default(groups))
if (missing(patterns)) { patterns <- rbind(rep(0,K),0:(K-1)); colnames(patterns) <- unique.default(groups) }
groupsr <- groups2int(groups,patterns)

if (length(table(groupsr))!=ncol(patterns)) stop('patterns must have the same number of columns as the number of distinct groups')
if (ncol(patterns)!=K) stop('patterns must have number of columns equal to the number of distinct elements in groups')
if (sum(is.na(patterns))>0) stop('patterns cannot have any NA values')
if (sum(is.nan(patterns))>0) stop('patterns cannot have any NaN values')
if (sum(is.infinite(patterns))>0) stop('patterns cannot have any Inf values')
for (i in 1:nrow(patterns)) { patterns[i,] <- as.integer(as.integer(as.factor(patterns[i,]))-1) }
class(patterns) <- 'gagahyp'
ngrouppat <- as.integer(apply(patterns,1,'max')+1)
nclust <- as.integer(nclust)
npat <- as.integer(nrow(patterns))

# Initialize parameters by method of moments
if (nclust==1) probclusini <- 1
if (missing(parini)) {
  if (trace) cat('Initializing parameters...')
  aest <- rowMeans(x)^2/((rowMeans(x^2)-rowMeans(x)^2)*ncol(x)/(ncol(x)-1)); lest <- 1/rowMeans(x)
  sel <- (aest<quantile(aest,probs=.99,na.rm=TRUE)) & (lest<quantile(lest,probs=.99,na.rm=TRUE))
  aest <- aest[sel]; lest <- lest[sel]
  balphaini <- as.double(mean(aest)^2/var(aest,na.rm=TRUE)); nualphaini <- as.double(mean(aest,na.rm=TRUE))
  if (nclust==1) {
    a0ini <- as.double(mean(lest)^2/var(lest,na.rm=TRUE)); nuini <- as.double(mean(lest,na.rm=TRUE))
    probclusini <- as.double(1)
  } else {
    clusini <- kmeans(x=lest,centers=nclust)$cluster
    a0ini <- as.double(tapply(lest,clusini,'mean')^2/tapply(lest,clusini,'var'))
    nuini <- as.double(tapply(lest,clusini,'mean'))
    probclusini <- as.double(table(clusini)/length(clusini))
  }
  probpatini <- rep(1/nrow(patterns),nrow(patterns))
  if (trace) cat(' Done.\n')
} else {
  if (is.null(parini$a0) | is.null(parini$nu) | is.null(parini$balpha) | is.null(parini$nualpha) | is.null(parini$probpat)) stop('some components of parini are empty')
  a0ini <- parini$a0; nuini <- parini$nu; balphaini <- parini$balpha; nualphaini <- parini$nualpha
  if (nclust==1) { probclusini <- 1 } else { if (is.null(parini$probclus)) { stop('component probclus of parini is empty') } else { probclusini <- parini$probclus/sum(parini$probclus) } }
  probpatini <- parini$probpat/sum(parini$probpat)
}


if (method=='EM' | method=='quickEM') {

  if (method=='quickEM') B <- 1
  balpha <- nualpha <- double(1); alpha0 <- nu <- probclus <- double(nclust); prob <- double(npat)
  lhood <- double(1); trace <- as.integer(trace)
  a0ini <- as.double(a0ini); nuini <- as.double(nuini); balphaini <- as.double(balphaini); nualphaini <- as.double(nualphaini)
  probclusini <- as.double(probclusini); probpatini <- as.double(probpatini)
  groupsr <- as.integer(groupsr); K <- as.integer(K); equalcv <- as.integer(equalcv)
  nclust <- as.integer(nclust); npat <- as.integer(npat); ngrouppat <- as.integer(ngrouppat)
  gapprox <- as.integer(gapprox); trace <- as.integer(trace)
  
  z <- .C("fitEM_ggC",alpha0=alpha0,nu=nu,balpha=balpha,nualpha=nualpha,probclus=probclus,prob=prob,lhood=lhood,as.integer(B),a0ini,nuini,balphaini,nualphaini,probclusini,probpatini,as.integer(nrow(x)),as.integer(ncol(x)),as.double(t(x)),groupsr,K,equalcv,nclust,npat,as.integer(t(patterns)),ngrouppat,gapprox,trace)
  
  parest <- c(alpha0=z$alpha0,nu=z$nu,balpha=z$balpha,nualpha=z$nualpha,probclus=z$probclus,probpat=z$prob)
  gg.fit <- list(parest=parest,mcmc=as.mcmc(NA),lhood=z$lhood,equalcv=equalcv,nclust=nclust,patterns=patterns,method=method)
  class(gg.fit) <- 'gagafit'
  return(gg.fit)

} else if (method=='MH' | method=='Gibbs' | method=='SA') {

  if (trace) cat('Refining initial estimates...')
  eps <- 1; i <- 1
  while ((eps>.001) && (i<=B)) {
    probnew <- colMeans(ppGG(x,groups,a0ini,nuini,balphaini,nualphaini,equalcv,probclusini,probpatini,patterns)$pp)
    eps <- max(abs(probnew-probpatini))
    probpatini <- probnew
    i <- i+1
  }
  if (trace) cat(' Done.\n')

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
  balpha <- nualpha <- double(B); alpha0 <- nu <- probclus <- double(B*nclust); prob <- double(B*npat)
  lhood <- double(B); trace <- as.integer(trace)

  if (method=='Gibbs') {
    z <- .C("fit_ggC",alpha0=alpha0,nu=nu,balpha=balpha,nualpha=nualpha,probclus=probclus,prob=prob,lhood=lhood,as.integer(B),as.double(a.alpha0),as.double(b.alpha0),as.double(a.nu),as.double(b.nu),as.double(a.balpha),as.double(b.balpha),as.double(a.nualpha),as.double(b.nualpha),as.double(p.probclus),as.double(p.probpat),a0ini,nuini,balphaini,nualphaini,probclusini,probpatini,as.integer(nrow(x)),as.integer(ncol(x)),as.double(t(x)),as.integer(groupsr),K,as.integer(equalcv),nclust,npat,as.integer(t(patterns)),ngrouppat,as.integer(gapprox),trace)
  } else if (method=='MH' | method=='SA') {
    acprop <- as.double(0)
    h.alpha0 <- h.nu <- h.balpha <- h.nualpha <- h.rho <- h.prob <- as.double(0)
    if (method=='MH') { cmethod <- as.integer(1); Bgibbs <- as.integer(50) } else { cmethod <- as.integer(0); Bgibbs <- as.integer(20) }
    z <- .C("fitMH_ggC",acprop=acprop,alpha0=alpha0,nu=nu,balpha=balpha,nualpha=nualpha,probclus=probclus,prob=prob,lhood=lhood,as.integer(B),as.double(a.alpha0),as.double(b.alpha0),as.double(a.nu),as.double(b.nu),as.double(a.balpha),as.double(b.balpha),as.double(a.nualpha),as.double(b.nualpha),as.double(p.probclus),as.double(p.probpat),a0ini,nuini,balphaini,nualphaini,probclusini,probpatini,as.integer(nrow(x)),as.integer(ncol(x)),as.double(t(x)),as.integer(groupsr),K,as.integer(equalcv),nclust,npat,as.integer(t(patterns)),ngrouppat,as.integer(gapprox),trace,cmethod,Bgibbs,h.alpha0,h.nu,h.balpha,h.nualpha,h.rho,h.prob)
  }
  if (trace) cat('Done.\n')

  gg.fit <- list(parest=NA,mcmc=as.mcmc(data.frame(alpha0=matrix(z$alpha0,nrow=B,ncol=nclust,byrow=TRUE),nu=matrix(z$nu,nrow=B,ncol=nclust,byrow=TRUE),balpha=z$balpha,nualpha=z$nualpha,probclus=matrix(z$probclus,nrow=B,ncol=nclust,byrow=TRUE),probpat=matrix(z$prob,nrow=B,byrow=TRUE))),lhood=z$lhood,equalcv=equalcv,nclust=nclust,patterns=patterns,method=method)
  class(gg.fit) <- 'gagafit'
  return(gg.fit)

}

}
