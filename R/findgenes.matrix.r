findgenes.matrix <- function(fit,x,groups,fdrmax=.05,parametric=TRUE,B=500) {
  cf <- as.double(2)
  if (ncol(fit)<2) stop('fit must be a matrix with >=2 columns. The 1st column must contain the post prob for the null hypothesis')
  nsel <- nrow(fit); sel <- as.integer((1:nsel)-1)
  util <- as.integer(1)
  u <- double(1); d <- integer(nrow(fit)); fdr <- fnr <- power <- double(1); threshold <- double(1)
  z <- .C("utgene_parC",u=u,d=d,fdr=fdr,fnr=fnr,power=power,threshold=threshold,util=util,as.double(cf),nsel,sel,as.double(t(fit)),as.integer(ncol(fit)),as.double(fdrmax))
  fdr <- fdrest <- fdrpar <- z$fdr
  return(list(truePos=z$u,d=z$d,fdr=fdr,fdrpar=fdrpar,fdrest=fdrest,fnr=z$fnr,power=z$power,threshold=z$threshold))
}
