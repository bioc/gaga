sim.gg <- function(n,m,p.de=.1,a0,nu,balpha,nualpha,probclus=1,a=NA,l=NA,useal=FALSE) {
  # Simulates data from the GaGa model with 2 groups
  # Input:
  # - n: number of genes
  # - m: number of observations per group
  # - p.de: proportion of genes differentially expressed
  # - a0, nu: prior for mean is IG(a0,a0/nu).
  # - balpha, nualpha: prior for alpha is G(balpha,balpha/nualpha)
  # Output:
  # - x: matrix with simulated observations for group 1 in first m columns and group 2 in last m columns
  # - l: matrix with lambda parameters used to generate x and y in column 1 and 2, respectively

  if (n<=0) stop("Number of genes must be positive")
  if (sum(m<0)>0) stop("Number of observations per group must be positive")

  if (useal==FALSE) {
    if (p.de<0 & p.de>1) stop("proportion of differentially expressed genes must be between 0 and 1")
    if (balpha<=0 | nualpha<=0) stop("balpha and nualpha must be >0")
    if (min(a0)<=0 | min(nu)<=0) stop(cat("a0 and nu must be >0 but a0=",a0,", nu=",nu,"was specified \n"))

    a0 <- rep(a0,round(probclus*n)); nu <- rep(nu,round(probclus*n))
    balpha <- rep(balpha,length(a0)); nualpha <- rep(nualpha,length(nu))
 
    a <- l <- matrix(NA,nrow=n,ncol=length(m))
    a[,1] <- rgamma(n,balpha,balpha/nualpha); l[,1] <- 1/rgamma(n,a0,a0/nu)
    if (ncol(a)>1) {
      l[,2] <- l[,1]; a[,2] <- a[,1]
      if ((round((1-p.de)*n)+1)<=n) {   #generate parameter values for DE genes
        sel <- sample(1:n,round(p.de*n),replace=FALSE)
        a[sel,2] <- rgamma(round(p.de*n),balpha[sel],balpha[sel]/nualpha[sel])
        l[sel,2] <- 1/rgamma(round(p.de*n),a0[sel],a0[sel]/nu[sel])
      }
    }
  }
  x <- matrix(rgamma(m[1]*n,a[,1],a[,1]/l[,1]),nrow=n,ncol=m[1])
  i <- 2
  while (i<=length(m)) {
    x <- cbind(x,matrix(rgamma(m[i]*n,a[,i],a[,i]/l[,i]),nrow=n,ncol=m[i]))
    i <- i+1
  }

  return(list(x=x,a=a,l=l))
}
