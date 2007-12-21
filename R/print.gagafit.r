print.gagafit <- function(x,...) {

if (x$nclust==1) cat("GaGa hierarchical model.") else cat("MiGaGa hierarchical model (",round(x$nclust)," clusters.",sep="")
if (x$method=='EBayes') cat(" Fit via empirical Bayes\n") else cat(" Fit via MCMC (",nrow(x$mcmc)," iterations kept)\n",sep="")

if (is.null(x$pp)) {
  cat("  ",ncol(x$patterns)," groups, ",nrow(x$patterns)," hypotheses (expression patterns)\n\n",sep="") } else {
  cat("  ",nrow(x$pp)," genes, ",ncol(x$patterns)," groups, ",nrow(x$patterns)," hypotheses (expression patterns)\n\n",sep="")
}
cat("The expression patterns are\n")
print(x$patterns)
cat("\n")
if (is.null(x$parest)==FALSE) {
  cat("Hyper-parameter estimates\n\n")
  par <- getpar(x)
  cat("  ",names(x$parest)[1:(2+2*x$nclust)],"\n")
  cat("  ",round(par$a0,3),round(par$nu,3),round(par$balpha,3),round(par$nualpha,3),"\n\n")
  cat("  ",names(x$parest)[(2+2*x$nclust+1):(2+3*x$nclust)],"\n")
  cat("  ",round(par$probclus,3),"\n\n")
  cat("  ",names(x$parest)[-1:-(2+3*x$nclust)],"\n")
  cat("  ",round(par$probpat,3),"\n")
}
if (is.null(x$pp)) { cat("Posterior probabilities not available. Run function parest.\n") }
if (is.null(x$parest)) { cat("Hyper-parameter estimates not available. Run function parest.\n") }
}
