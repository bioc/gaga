getpar.nnfit <- function(fit, ...) { fit$parest }

print.nnfit <- function(x,...) {
cat("Normal-Normal hierarchical model.")
cat(" Fit via Expectation-Maximization\n")
cat("Assumed different variance for each gene (NNGV model in EBarrays) \n")
#
if (is.null(x$pp)) {
  cat("  ",ncol(x$patterns)," groups, ",nrow(x$patterns)," hypotheses (expression patterns)\n\n",sep="") } else {
  cat("  ",nrow(x$pp)," genes, ",ncol(x$patterns)," groups, ",nrow(x$patterns)," hypotheses (expression patterns)\n\n",sep="")
}
cat("The expression patterns are\n")
probpat <- getpar(x)
print(x$patterns,probpat[6:7])
cat("\n")
cat("Hyper-parameter estimates\n\n")
par <- getpar(x)
cat(paste(paste(names(par[1:4]),round(par[1:4],3),sep='='),collapse=' '),'\n')
}
