simGG <- function(n,m,p.de=.1,a0,nu,balpha,nualpha,equalcv=TRUE,probclus=1,a=NA,l=NA,useal=FALSE) {
  # Simulates data from the GaGa model with several groups
  if (n<=0) stop("Number of genes must be positive")
  if (sum(m<0)>0) stop("Number of observations per group must be positive")
  if (!missing(a)) { if (length(m)!=ncol(a)) stop("length(m) must be equal to ncol(a)") }
  if (!missing(l)) { if (length(m)!=ncol(l)) stop("length(m) must be equal to ncol(l)") }
  
  if (useal==FALSE) {
    if (p.de<0 & p.de>1) stop("proportion of differentially expressed genes must be between 0 and 1")
    if (balpha<=0 | nualpha<=0) stop("balpha and nualpha must be >0")
    if (min(a0)<=0 | min(nu)<=0) stop(cat("a0 and nu must be >0 but a0=",a0,", nu=",nu,"was specified \n"))

    a0 <- rep(a0,round(probclus*n)); nu <- rep(nu,round(probclus*n))
    balpha <- rep(balpha,length(a0)); nualpha <- rep(nualpha,length(nu))
 
    a <- l <- matrix(NA,nrow=n,ncol=length(m))
    a[,1] <- rgamma(n,balpha,balpha/nualpha); l[,1] <- 1/rgamma(n,a0,a0/nu)
    if (ncol(a)>1) {
      if ((round((1-p.de)*n)+1)<=n) {   #generate parameter values for DE genes
        sel <- sample(1:n,round(p.de*n),replace=FALSE)
        for (i in 2:length(m)) {
          l[,i] <- l[,1]; a[,i] <- a[,1]
          if (!equalcv) a[sel,i] <- rgamma(round(p.de*n),balpha[sel],balpha[sel]/nualpha[sel])
          l[sel,i] <- 1/rgamma(round(p.de*n),a0[sel],a0[sel]/nu[sel])
        }
      }
    }
  }
  x <- matrix(rgamma(m[1]*n,a[,1],a[,1]/l[,1]),nrow=n,ncol=m[1])
  i <- 2
  while (i<=length(m)) {
    x <- cbind(x,matrix(rgamma(m[i]*n,a[,i],a[,i]/l[,i]),nrow=n,ncol=m[i]))
    i <- i+1
  }

  metadata <- data.frame(labelDescription='Group that each array belongs to',row.names='group')
  group <- paste('group',1:length(m))
  pheno <- new("AnnotatedDataFrame", data=data.frame(group=rep(group,m)), dimLabels=c("rowNames", "columnNames"), varMetadata=metadata)
  sampleNames(pheno) <- paste("Array",1:nrow(pheno))
  metadata <- data.frame(labelDescription=c(paste('alpha parameter for array',1:length(m)),paste('mean parameter for array',1:length(m))),row.names=c(paste('alpha',1:length(m),sep='.'),paste('mean',1:length(m),sep='.')))
  fdata <- new("AnnotatedDataFrame", data=data.frame(a,l),varMetadata=metadata)
  sampleNames(fdata) <- paste("Gene ",1:n)
  varLabels(fdata) <- c(paste('alpha',1:length(m),sep='.'),paste('mean',1:length(m),sep='.'))
  experimentData <- new("MIAME", title = "Dataset simulated with simGG", abstract = "This dataset contains simulated expression data from a GaGa model. The expression data can be accessed via exprs(object), and the parameter values used to generate the data through fData(object)")
  x <- new("ExpressionSet", phenoData=pheno, featureData=fdata, exprs = x, experimentData = experimentData)
  
  return(x)
}
