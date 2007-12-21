fitfreq.gg <- function(x,groups,patterns,nclust=1,a0ini,nuini,balphaini,nualphaini,probpatini,iter.max=100,trace=FALSE) {
# Fits Generalized Gamma/Gamma model by maximum likelihood (uses nlminb)
# Input:
# - x: matrix with gene expression measurements for all groups
# - groups: vector of length ncol(x) indicating what group does each column in x correspond to
# - patterns: matrix indicating which groups are put together under each pattern, i.e. the hypotheses to consider for each gene. Defaults to two hypotheses: null hypothesis of all groups being equal and full alternative of all groups being different
# - a0ini,nuini,balphaini,nualphaini,probpatini: initial values for hyper-parameters in optimization
# - iter.max: maximum number of iterations in the optimization process
# - trace: if trace==TRUE the optimization progress is displayed every iteration
# Output: a list with the following components
# - a0, nu, balpha, nualpha, prob: hyper-parameter estimates
# - pp: posterior probabilities of each expression pattern

gapprox <- TRUE
cluslist <- c(0,-1); probclus <- 1 #EM only implemented for 1 cluster

ilogit <- function(x) { return(1/(1+exp(-x))) }
logl <- function(th) {
a0 <- as.double(exp(th[1:nclust])); nu <- as.double(exp(th[(nclust+1):(2*nclust)]))
balpha <- as.double(exp(th[2*nclust+1])); nualpha <- as.double(exp(th[2*nclust+2]))
z <- pp.gg(x,groups,a0,nu,balpha,nualpha,probclus,prob,patterns)
return(-z$lhood)
}


#EM algorithm
a0 <- a0ini; nu <- nuini; balpha <- balphaini; nualpha <- nualphaini; prob <- probpatini
th <- double(2*nclust+2); th[1:nclust] <- log(a0); th[(nclust+1):(2*nclust)] <- log(nu)
th[2*nclust+1] <- log(balpha); th[2*nclust+2] <- log(nualpha)
i <- 1; thnorm <- lnorm <- lcur <- ldif <- 1
if (trace) cat("Iter",format('a0',width=8,justify='right'),format('nu',width=8,justify='right'),format('balpha',width=8,justify='right'),format('nualpha',width=8,justify='right'),format('prob',width=length(prob)*8+1,justify='right'),'log-likelihood',"\n")
while ((i<=iter.max) & ((thnorm>1e-3) | (lnorm>1e-5) | (ldif>1e-5))) {
z <- pp.gg(x,groups,a0,nu,balpha,nualpha,probclus,prob,patterns)
probold <- prob; prob <- colMeans(z$pp)
mstep <- nlminb(start=th,objective=logl,control=list(iter.max=50,trace=0,abs.tol=1e-5,rel.tol=1e-5,x.tol=1e-5))
thnorm <- sqrt( sum((mstep$par-th)^2) + sum((ilogit(prob[-1])-ilogit(probold[-1]))^2) )
lnorm <- abs(mstep$objective-lcur)/abs(lcur); ldif <- abs(mstep$objective-lcur)
th <- mstep$par; lcur <- mstep$objective
a0 <- exp(th[1:nclust]); nu <- exp(th[(nclust+1):(2*nclust)])
balpha <- exp(th[2*nclust+1]); nualpha <- exp(th[2*nclust+2])
if (trace>0) { cat(format(i,width=4,justify='right'),round(a0,6),round(nu,6),round(balpha,6),round(nualpha,6),format(round(prob,6),width=8),format(-lcur,nsmall=3),"\n") }
i <- i+1
}

z <- pp.gg(x,groups,a0,nu,balpha,nualpha,probclus,prob,patterns)
return(list(alpha0=a0,nu=nu,balpha=balpha,nualpha=nualpha,probpat=prob,lhood=z$lhood))
}
