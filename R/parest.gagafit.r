parest.gagafit <- function(gg.fit,x,groups,burnin,alpha=.05) {
# Parameter estimates and posterior probabilities of differential expression for GaGa and MiGaGa model
# Input:
# - gg.fit: GaGa model fit as returned by fit.gg
# - x:  matrix with gene expression measurements for all groups
# - groups: vector of length ncol(x) indicating what group does each column in x correspond to
# - burnin: if gg.fit was fit via MCMC, burnin indicates number of MCMC iterations to be discarded
# - alpha: credibility interval with 1-alpha posterior probability is returned
# Output:
# - a0, nu, balpha, nualpha, probclus, probpat: hyper-parameter estimates (posterior mean for Bayes, maximum likelihood estimate for EBayes)
# - ci: posterior credibility intervals for the hyper-parameter estimates (only available for Bayes fit)
# - nclust: number of clusters
# - pp: posterior probabilities of each expression pattern
# - patterns: matrix indicating which groups are put together under each pattern, as passed to fit.gg
# - dic: DIC (only returned for Bayes fit)

if (missing(x)) stop('argument x must be specified')
if (missing(groups)) stop('argument groups must be specified')
if (ncol(x)!=length(groups)) stop('length(groups) must be equal to the number of columns in x')

nclust <- gg.fit$nclust
if (gg.fit$method=='EBayes') {
  a0 <- gg.fit$parest[1]; nu <- gg.fit$parest[2]
  balpha <- gg.fit$parest[3]; nualpha <- gg.fit$parest[4]
  probclus <- 1; probpat <- gg.fit$parest[6:7]
  ci<-list(a0=NA,nu=NA,balpha=NA,nualpha=NA,probclus=NA,probpat=NA)
  pp <- pp.gg(x,groups,a0=a0,nu=nu,balpha=balpha,nualpha=nualpha,probclus=probclus,probpat=probpat,patterns=gg.fit$patterns)
  dic <- NA
} else {
  if (missing(burnin)) {
    warning('burnin not specified, discarding 25% of the MCMC samples')
    gg.fit$mcmc <- gg.fit$mcmc[-1:-(.25*nrow(gg.fit$mcmc)),]
    gg.fit$lhood <- gg.fit$lhood[-1:-(.25*nrow(gg.fit$mcmc))]
    lhood <- mean(gg.fit$lhood)
  } else {
    if (burnin>=nrow(gg.fit$mcmc)) stop('burnin must be smaller than the number of MCMC samples')
    gg.fit$mcmc <- gg.fit$mcmc[-1:-burnin,]
    gg.fit$lhood <- gg.fit$lhood[-1:-burnin]
    lhood <- mean(gg.fit$lhood)
  }

  if (nclust>1) {
    a0 <- colMeans(gg.fit$mcmc[,1:nclust])
    nu <- colMeans(gg.fit$mcmc[,(nclust+1):(2*nclust)])
    balpha <- mean(gg.fit$mcmc[,2*nclust+1])
    nualpha <- mean(gg.fit$mcmc[,2*nclust+2])
    probclus <- colMeans(gg.fit$mcmc[,(2*nclust+3):(3*nclust+2)])
    probpat <- colMeans(gg.fit$mcmc[,-1:-(3*nclust+2)])
    a0.ci <- apply(gg.fit$mcmc[,1:nclust],2,quantile,probs=c(alpha/2,1-alpha/2))
    nu.ci <- apply(gg.fit$mcmc[,(nclust+1):(2*nclust)],2,quantile,probs=c(alpha/2,1-alpha/2))
    balpha.ci <- quantile(gg.fit$mcmc[,2*nclust+1],probs=c(alpha/2,1-alpha/2))
    nualpha.ci <- quantile(gg.fit$mcmc[,2*nclust+2],probs=c(alpha/2,1-alpha/2))
    probclus.ci <- apply(gg.fit$mcmc[,(2*nclust+3):(3*nclust+2)],2,quantile,probs=c(alpha/2,1-alpha/2))
    probpat.ci <- apply(gg.fit$mcmc[,-1:-(3*nclust+2)],2,quantile,probs=c(alpha/2,1-alpha/2))
  } else {
    a0 <- mean(gg.fit$mcmc[,1:nclust])
    nu <- mean(gg.fit$mcmc[,(nclust+1):(2*nclust)])
    balpha <- mean(gg.fit$mcmc[,(2*nclust+1):(3*nclust)])
    nualpha <- mean(gg.fit$mcmc[,(3*nclust+1):(4*nclust)])
    probclus <- mean(gg.fit$mcmc[,(4*nclust+1):(5*nclust)])
    probpat <- colMeans(gg.fit$mcmc[,-1:-(5*nclust)])
    a0.ci <- quantile(gg.fit$mcmc[,1:nclust],probs=c(alpha/2,1-alpha/2))
    nu.ci <- quantile(gg.fit$mcmc[,(nclust+1):(2*nclust)],probs=c(alpha/2,1-alpha/2))
    balpha.ci <- quantile(gg.fit$mcmc[,2*nclust+1],probs=c(alpha/2,1-alpha/2))
    nualpha.ci <- quantile(gg.fit$mcmc[,2*nclust+2],probs=c(alpha/2,1-alpha/2))
    probclus.ci <- c(1,1)
    probpat.ci <- apply(gg.fit$mcmc[,-1:-(3*nclust+2)],2,quantile,probs=c(alpha/2,1-alpha/2))
  }

  gg.fit$parest <- c(a0=a0,nu=nu,balpha=balpha,nualpha=nualpha,probclus=probclus,probpat)
  ci<-list(a0=a0.ci,nu=nu.ci,balpha=balpha.ci,nualpha=nualpha.ci,probclus=probclus.ci,probpat=probpat.ci)
  pp <- pp.gg(x,groups,a0=a0,nu=nu,balpha=balpha,nualpha=nualpha,probclus=probclus,probpat=probpat,patterns=gg.fit$patterns)
  dic <- -2*(2*lhood-pp$lhood)
}

gg.fit$ci <- ci; gg.fit$pp <- pp$pp; gg.fit$dic <- dic
return(gg.fit)
}
