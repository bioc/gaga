forwsimDiffExpr <- function(gg.fit,x,groups,ngenes,maxBatch,batchSize,fdrmax=.05,genelimit,v0thre=1,B=100,Bsummary=100,trace=TRUE,randomSeed) {

if (missing(batchSize)) stop('batchSize must be specified')
batchSize <- as.integer(batchSize)
if (missing(maxBatch)) stop('maxBatch must be specified')
if (round(maxBatch)!=maxBatch | maxBatch<1) stop('maxBatch must be an integer >=1')
if (round(B)!=B | B<1) stop('B must be an integer >=1')

gapprox <- TRUE
patterns <- gg.fit$patterns
if (missing(x)) {
  simprior <- as.integer(1)
  x <- double(0); groupsr <- integer(0)
  K <- as.integer(nrow(patterns))
  nrowx <- as.integer(ngenes)
  ncolx <- as.integer(0)
  if (missing(genelimit)) { genelimit <- ngenes }
} else {
  simprior <- as.integer(0)
  if (is(x, "exprSet") | is(x,"ExpressionSet")) {
    if (is.character(groups) && length(groups)==1) { groups <- as.factor(pData(x)[, groups]) }
    x <- exprs(x)
  } else if (!is(x,"data.frame") & !is(x,"matrix")) { stop("x must be an exprSet, data.frame or matrix") }
  if (ncol(x)!=length(groups)) stop('Argument groups must be of length equal to the number of columns in x')
  if (missing(genelimit)) { genelimit <- nrow(x); }
  groupsr <- groups2int(groups,patterns)
  K <- as.integer(max(groupsr)+1)
  if (ncol(patterns)!=K) stop('patterns must have number of columns equal to the number of distinct elements in groups')
  nrowx <- as.integer(nrow(x))
  ncolx <- as.integer(ncol(x))
}

#Set up some C routine variables that are not currently needed
stepsum <- 1 #summary statistic is computed up to 'stepsum' steps ahead in time
cf <- double(2); cf.sam <- as.double(0) #utility coeff. These values will make the C routine to ignore them.
Jini <- 0 #initial sample size

B <- as.integer(B); Bsummary <- as.integer(Bsummary); maxBatch <- as.integer(maxBatch); Jini <- as.integer(Jini)
util <- as.integer(1)
for (i in 1:nrow(patterns)) { patterns[i,] <- as.integer(as.integer(as.factor(patterns[i,]))-1) }
ngrouppat <- as.integer(apply(patterns,1,'max')+1)
par <- getpar(gg.fit)
a0 <- as.double(par$a0); nu <- as.double(par$nu)
balpha <- as.double(par$balpha); nualpha <- as.double(par$nualpha)
probclus <- as.double(par$probclus); probpat <- as.double(par$probpat)
fdrmax <- as.double(fdrmax)
simid <- j <- integer(B*(maxBatch-Jini)); u <- fdr <- fnr <- power <- double(B*(maxBatch-Jini))
summary <- double(stepsum*B*(maxBatch-Jini))

#Set random number generator seed
if (missing(randomSeed)) { randomSeed <- as.numeric(Sys.time()) }
if (randomSeed<=0) stop('randomSeed must be >0')
randomSeed <- randomSeed + 2 #must be >2 to avoid problems with C random number seed
randomSeed <- as.integer(c(randomSeed,as.integer(randomSeed/2)) %% 10^6)

cat("Starting forward simulation...\n")
z <- .C("forwsim_geneC",simid=simid,j=j,u=u,fdr=fdr,fnr=fnr,power=power,summary=summary,B,Bsummary,as.integer(stepsum),maxBatch,Jini,batchSize,util,cf,cf.sam,as.integer(genelimit),as.double(v0thre),simprior,nrowx,ncolx,as.double(t(x)),groupsr,K,a0,nu,balpha,nualpha,as.integer(gg.fit$equalcv),as.integer(length(probclus)),probclus,probpat,as.integer(nrow(patterns)),as.integer(t(patterns)),ngrouppat,fdrmax,as.integer(trace),as.integer(gapprox),randomSeed)

if (stepsum) { summary <- matrix(z$summary,ncol=stepsum,byrow=TRUE) } else { summary <- rep(NA,length(z$simid)) }
ans <- data.frame(simid=z$simid,time=z$j,u=round(z$u,5),fdr=round(z$fdr,5),fnr=round(z$fnr,5),power=round(z$power,5),summary=round(summary,5))

if (simprior) {
  u0 <- list(truePos=0,fdr=1,fnr=1,power=0)
} else {
  u0 <- findgenes(gg.fit=gg.fit,x=x,groups=groups,fdrmax=fdrmax,parametric=TRUE)
}
s0 <- mean(ans$u[ans$time==min(ans$time)])-u0$truePos
ans0 <- data.frame(simid=unique(ans$simid),time=min(ans$time)-1,u=u0$truePos,fdr=u0$fdr,fnr=u0$fnr,power=u0$power,summary=s0)
ans <- rbind(ans0,ans); ans <- ans[order(ans$simid,ans$time),]
ans$summary[ans$summary==-9999] <- NA

return(ans)
}
