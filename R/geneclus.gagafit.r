geneclus.gagafit <- function(gg.fit,method='posprob') {

  d <- integer(nrow(gg.fit$pp)); ppat <- double(nrow(gg.fit$pp))
  nsel <- as.integer(nrow(gg.fit$pp)); sel <- as.integer((1:nsel)-1)
  npat <- as.integer(nrow(gg.fit$patterns))
  probpat <- as.double(getpar(gg.fit)$probpat)
  if (method=='likelihood') { gg.fit$pp <- t(t(gg.fit$pp)/probpat); gg.fit$pp <- gg.fit$pp/rowSums(gg.fit$pp) }
  v <- as.double(t(gg.fit$pp))

  z <- .C("geneclus",d=d,ppat=ppat,nsel,sel,v,npat)
  z <- list(d=z$d,posprob=z$ppat,patterns=gg.fit$patterns)
  class(z) <- 'gagaclus'
  return(z)

}
