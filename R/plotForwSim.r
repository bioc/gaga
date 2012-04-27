plotForwSim <- function(fs,xlab="Number of additional batches per group",ylab="Expected increase in True Positives",col='gray',lty=1,...) {
  sel <- fs$simid==min(fs$simid)
  plot(fs$time[sel],fs$summary[sel],xlab=xlab,ylab=ylab,type='l',col=col,lty=lty,...)
  for (i in (min(fs$simid)+1):max(fs$simid)) {
    sel <- fs$simid==i
    lines(fs$time[sel],fs$summary[sel],col=col,lty=lty)
  }
}
