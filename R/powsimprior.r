powsimprior <- function(fit, m, ngenes, fdrmax=.05, B=1000, mc.cores=1) { UseMethod("powsimprior") }

powsimprior.nnfit <- function(fit, m, ngenes, fdrmax=.05, B=1000, mc.cores=1) {
  est <- getpar(fit)
  probpat <- est[grep('probpat',names(est))]
  if (length(probpat)>2) stop('Only two patterns currently implemented')
  mu0 <- est['mu0']; tau0 <- sqrt(est['tau02']); v0 <- est['v0']; sigma0 <- sqrt(est['sigma02'])
  fit$pp <- matrix(rep(probpat,each=nrow(fit$pp)),nrow=nrow(fit$pp),ncol=length(probpat))
  f <- function(simid) {
    xnew <- simNN(n=ngenes, m=m, p.de=probpat['probpat2'],mu0=mu0,tau0=tau0,v0=v0,sigma0=sigma0)
    fitnew <- updateNNfit(fit=fit,xnew=xnew,probpat=probpat)
    findgenes(fitnew,fdrmax=fdrmax)$truePos
  }
  if (mc.cores>1) {
    if ('multicore' %in% loadedNamespaces()) {
      ans <- multicore::mclapply(1:B, f, mc.cores=mc.cores, mc.preschedule=TRUE)
    } else stop('multicore library has not been loaded!')
  } else {
    ans <- lapply(1:B, f)
  }
  ans <- unlist(ans)
  return(list(m=mean(ans), s=sd(ans)/sqrt(B)))
}
