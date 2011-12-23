simLNN <- function(n, m, p.de=0.1, mu0, tau0, v0, sigma0) {
  # Simulate x_ij ~ logN(mu_ij,sigma_j); mu_ij ~ N(mu0,tau0^2); sigma_j^2 ~ IG(.5*v0, .5*v0*sigma0^2)
  ans <- simNN(n=n, m=m, p.de=p.de, mu0=mu0, tau0=tau0, v0=v0, sigma0=sigma0)
  exprs(ans) <- exp(exprs(ans))
  return(ans)
}

simNN <- function(n, m, p.de=0.1, mu0, tau0, v0, sigma0) {
  # Simulate x_ij ~ N(mu_ij,sigma_j); mu_ij ~ N(mu0,tau0^2); sigma_j^2 ~ IG(.5*v0, .5*v0*sigma0^2)
  if (n <= 0) stop("Number of genes must be positive")
  if (sum(m < 0) > 0) stop("Number of observations per group must be positive")
  if (p.de < 0 & p.de > 1) stop("proportion of differentially expressed genes must be between 0 and 1")
  if (tau0 <= 0) stop("tau0 must be >0")
  if (v0 <= 0) stop("v0 must be >0")
  if (sigma0 <= 0) stop("sigma0 must be >0")
  #Simulate parameter values
  d <- runif(n)<p.de
  mu <- matrix(NA,nrow=n,ncol=length(m))
  mu[,1] <- rnorm(n,mu0,sd=tau0)
  mu[d==0,] <- mu[d==0,1]
  sel <- d==1; nsel <- sum(sel); if (ncol(mu)>1) for (j in 2:ncol(mu)) mu[sel,j] <- rnorm(nsel,mu0,sd=tau0)
  sigma <- sqrt(1/rgamma(n,.5*v0, .5*v0*sigma0^2))
  #Simulate observations
  x <- matrix(NA,nrow=n,ncol=sum(m))
  stcol <- c(1,1+cumsum(m))
  for (i in 1:length(m)) x[,stcol[i]:(stcol[i+1]-1)] <- rnorm(n*m[i],mu[,i],sd=sigma)
  #Format as ExpressionSet
  fdata <- data.frame(mu,sigma); names(fdata)[1:length(m)] <- paste('mu',1:length(m),sep='')
  fdata <- new("AnnotatedDataFrame",fdata)
  pdata <- data.frame(group=rep(paste('group',1:length(m),sep=''),m))
  pdata <- new("AnnotatedDataFrame",pdata)
  ans <- new("ExpressionSet",exprs=x,phenoData=pdata,featureData=fdata)
  return(ans)
}
