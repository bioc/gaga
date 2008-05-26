print.gagafit <- function(x,...) {

if (x$nclust==1) cat("GaGa hierarchical model.") else cat("MiGaGa hierarchical model (",round(x$nclust)," clusters.",sep="")
if (x$method=='EM') cat(" Fit via Expectation-Maximization\n") else if (x$method=='quickEM') cat("Fit via quick Expectation-Maximization\n") else if (x$method=='Gibbs') cat(" Fit via Gibbs sampling (",nrow(x$mcmc)," iterations kept)\n",sep="") else if (x$method=='MH') cat(" Fit via Metropolis-Hastings sampling (",nrow(x$mcmc)," iterations kept)\n",sep="") else if (x$method=='SA') cat(" Fit via Simulated Annealing (",nrow(x$mcmc)," iterations)\n",sep="")
if (x$equalcv) cat("Assumed constant CV across groups \n") else cat("Assumed varying CV across groups\n")

if (is.null(x$pp)) {
  cat("  ",ncol(x$patterns)," groups, ",nrow(x$patterns)," hypotheses (expression patterns)\n\n",sep="") } else {
  cat("  ",nrow(x$pp)," genes, ",ncol(x$patterns)," groups, ",nrow(x$patterns)," hypotheses (expression patterns)\n\n",sep="")
}
cat("The expression patterns are\n")
if (sum(is.na(x$parest))==0) { probpat <- getpar(x)$probpat } else { probpat <- NA }
print(x$patterns,probpat)
cat("\n")
if (!is.na(x$parest[1])) {
  cat("Hyper-parameter estimates\n\n")
  par <- getpar(x)
  cat("  ",names(x$parest)[1:(2+2*x$nclust)],"\n")
  cat("  ",round(par$a0,3),round(par$nu,3),round(par$balpha,3),round(par$nualpha,3),"\n\n")
  cat("  ",names(x$parest)[(2+2*x$nclust+1):(2+3*x$nclust)],"\n")
  cat("  ",round(par$probclus,3),"\n\n")
} else {
  cat("Hyper-parameter estimates not computed yet. Run function parest.\n")
}
if (is.null(x$pp)) { cat("Posterior probabilities not computed yet. Run function parest.\n") }
}
