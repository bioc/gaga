simnewsamples <- function(fit,groupsnew,sel,x,groups) { UseMethod("simnewsamples") }

simnewsamples.nnfit <- function(fit,groupsnew,sel,x,groups) {
# Simulate parameters from the posterior of NNGV, and new observations from the posterior predictive
if (is(x, "exprSet") | is(x,"ExpressionSet")) {
  if (is.character(groups) && length(groups)==1) { groups <- as.factor(pData(x)[, groups]) }
  x <- exprs(x)
} else if (!is(x,"data.frame") & !is(x,"matrix")) { stop("x must be an exprSet, data.frame or matrix") } 
patterns <- fit$patterns
groupsr <- groups2int(groups,patterns); K <- as.integer(max(groupsr)+1)
groupsnewr <- groups2int(groupsnew,patterns)
v <- fit$pp

# Checks
if ((max(groupsnewr)>max(groupsr)) | (min(groupsnewr)<min(groupsr))) stop('Groups indicated in groupsnew do not match with those indicated in groups')
if (ncol(x)!=length(groups)) stop('length(groups) must be equal to the number of columns in x')
if (!is.matrix(v)) stop('Argument v must be a matrix')
if (ncol(v)!=nrow(patterns)) stop('Argument v must have as many columns as rows has patterns')
if (nrow(x)!=nrow(v)) stop('Arguments x and v must have the same number of rows')
if (ncol(patterns)!=K) stop('patterns must have number of columns equal to the number of distinct elements in groups')
for (i in 1:nrow(patterns)) { patterns[i,] <- as.integer(as.integer(as.factor(patterns[i,]))-1) }
par <- getpar(fit)

sigmanew <- double(nrow(x))
munew <- matrix(NA,nrow=nrow(x),ncol=length(unique(groupsnewr))); colnames(munew) <- colnames(patterns)
xnew <- matrix(NA,nrow=nrow(x),ncol=length(groupsnewr))

#Draw delta
if (nrow(patterns)>1) {
  for (j in 2:ncol(v)) v[,j] <- v[,j]+v[,j-1]
  u <- runif(nrow(x)); d <- apply(u<v,1,function(z) which(z)[1])
} else {
  d <- rep(1,nrow(x))
}

#Draw sigma
a.sigma <- .5*(par['v0']+ncol(x))
b.sigma <- rep(.5*par['v0']*par['sigma02'], nrow(x))
simpat <- unique(d)
for (k in 1:length(simpat)) {
  rowsel <- d==k
  curgroups <- unique(patterns[k,])
  for (j in 1:length(curgroups)) {
    groupids <- colnames(patterns)[patterns[k,]==curgroups[j]]
    colsel <- groups %in% groupids
    xbar <- rowMeans(x[rowsel,colsel])
    b.sigma[rowsel] <- b.sigma[rowsel] + .5*rowSums((x[rowsel,colsel]-xbar)^2)
  }
}
sigmanew <- 1/sqrt(rgamma(nrow(x),a.sigma,b.sigma))

#Draw mu
for (k in 1:length(simpat)) {
  rowsel <- d==k
  curgroups <- unique(patterns[k,])
  for (j in 1:length(curgroups)) {
    groupids <- colnames(patterns)[patterns[k,]==curgroups[j]]
    colsel <- groups %in% groupids; ncolsel <- sum(colsel)
    xbar <- rowMeans(x[rowsel,colsel])
    v <- 1/(ncolsel/sigmanew[rowsel]^2 + 1/par['tau02'])
    m <- (ncolsel*xbar/sigmanew[rowsel]^2 + par['mu0']/par['tau02']) * v
    munew[rowsel, groupids] <- rnorm(sum(rowsel),m,sd=sqrt(v))
  }
}

#Draw x
design <- as.matrix(model.matrix(~ -1 + factor(groupsnewr)))
xnew <- t(design %*% t(munew)) + matrix(rnorm(nrow(x)*length(groupsnewr),0,sd=sigmanew),nrow=nrow(xnew))

#Create ExpressionSet
metadata <- data.frame(labelDescription='Group that each array belongs to',row.names='group')
pheno <- new("AnnotatedDataFrame", data=data.frame(group=groupsnew), dimLabels=c("rowNames", "columnNames"), varMetadata=metadata)
sampleNames(pheno) <- paste("Array",1:nrow(pheno))
#
metadata <- data.frame(labelDescription=c('Expression pattern','Standard deviation',paste('mean',colnames(patterns))),row.names=c('d','sigma',paste('mean',1:ncol(patterns),sep='.')))
fdata <- data.frame(d,sigmanew,munew)
fdata <- new("AnnotatedDataFrame", data=fdata,varMetadata=metadata)
sampleNames(fdata) <- paste("Gene",1:nrow(fdata))
varLabels(fdata) <- c('d','sigma',paste('mean',1:ncol(patterns),sep='.'))
#
experimentData <- new("MIAME", title = "Dataset simulated with simnewsamples", abstract = "This dataset contains expression data simulated from the posterior predictive distribution of a normal-normal model with generalized variances (see emfit in package EBarrays). The expression data can be accessed via exprs(object), and the parameter values used to generate the data through fData(object)")
#ans <- data.frame(xnew)
colnames(xnew) <- paste('Array',1:nrow(pheno)); rownames(xnew) <- paste('Gene',1:nrow(fdata))
xnew <- new("ExpressionSet", phenoData=pheno, featureData=fdata, exprs = xnew, experimentData = experimentData)
return(xnew)
}


simnewsamples.gagafit <- function(fit,groupsnew,sel,x,groups) {
# Simulate parameters from the posterior of a GaGa or MiGaGa model, and new observations from the posterior predictive

gapprox <- TRUE
if (is(x, "exprSet") | is(x,"ExpressionSet")) {
  if (is.character(groups) && length(groups)==1) { groups <- as.factor(pData(x)[, groups]) }
  x <- exprs(x)
} else if (!is(x,"data.frame") & !is(x,"matrix")) { stop("x must be an exprSet, data.frame or matrix") } 

patterns <- fit$patterns
v <- fit$pp

groupsr <- groups2int(groups,patterns); K <- as.integer(max(groupsr)+1)
groupsnewr <- groups2int(groupsnew,patterns)
#groupsnewr <- integer(length(groupsnew))
#for (i in 1:ncol(patterns)) { groupsnewr[groupsnew == colnames(patterns)[i]] <- i - 1 }
#groupsnewr <- as.integer(groupsnewr)

if ((max(groupsnewr)>max(groupsr)) | (min(groupsnewr)<min(groupsr))) stop('Groups indicated in groupsnew do not match with those indicated in groups')
if (ncol(x)!=length(groups)) stop('length(groups) must be equal to the number of columns in x')
if (!is.matrix(v)) stop('Argument v must be a matrix')
if (ncol(v)!=nrow(patterns)) stop('Argument v must have as many columns as rows has patterns')
if (nrow(x)!=nrow(v)) stop('Arguments x and v must have the same number of rows')

if (ncol(patterns)!=K) stop('patterns must have number of columns equal to the number of distinct elements in groups')
for (i in 1:nrow(patterns)) { patterns[i,] <- as.integer(as.integer(as.factor(patterns[i,]))-1) }
ngrouppat <- as.integer(apply(patterns,1,'max')+1)
par <- getpar(fit)
alpha0 <- as.double(par$a0); nu <- as.double(par$nu); balpha <- as.double(par$balpha)
nualpha <- as.double(par$nualpha)
nclust <- as.integer(fit$nclust); rho <- as.double(par$probclus)
sumx <- double(nrow(x)*sum(ngrouppat)); nobsx <- double(sum(ngrouppat))
if (fit$equalcv) {
  prodx <- double(nrow(x))
} else {
  prodx <- double(nrow(x)*sum(ngrouppat))
}
gapprox <- as.integer(gapprox)

if (missing(sel)) sel <- 1:nrow(x)
if (is.logical(sel)) sel <- (1:nrow(x))[sel]
sel <- as.integer(sel-1) #in C indices start at 0
nsel <- length(sel); nsamples <- length(groupsnewr)
xnew <- anew <- lnew <- double(nsel*nsamples); dnew <- integer(nsel*nsamples)

z <- .C("compute_sumxC",sumx=t(sumx),prodx=t(prodx),nobsx=nobsx,as.integer(fit$equalcv),nsel,sel,as.integer(sum(ngrouppat)),ncol(x),as.double(t(x)),groupsr,K,as.integer(nrow(patterns)),as.integer(t(patterns)),ngrouppat,as.integer(1))
sumx <- matrix(z$sumx,nrow=nrow(x),byrow=TRUE); prodx <- matrix(z$prodx,nrow=nrow(x),byrow=TRUE)
nobsx <- z$nobsx

z <- .C("simnewsamples_ggC",xnew=xnew,dnew=dnew,anew=anew,lnew=lnew,nsamples,groupsnewr,K,nsel,sel,alpha0,nu,balpha,nualpha,as.integer(fit$equalcv),nclust,rho,as.double(t(v)),as.integer(nrow(patterns)),as.integer(t(patterns)),ngrouppat,t(sumx),t(prodx),nobsx,gapprox)

metadata <- data.frame(labelDescription='Group that each array belongs to',row.names='group')
pheno <- new("AnnotatedDataFrame", data=data.frame(group=groupsnew), dimLabels=c("rowNames", "columnNames"), varMetadata=metadata)
sampleNames(pheno) <- paste("Array",1:nrow(pheno))

metadata <- data.frame(labelDescription=c(paste('Expression patterns for array',1:length(groupsnewr)),paste('alpha parameters for array',1:length(groupsnewr)),paste('mean parameters for array',1:length(groupsnewr))),row.names=c(paste('d',1:length(groupsnewr),sep='.'),paste('alpha',1:length(groupsnewr),sep='.'),paste('mean',1:length(groupsnewr),sep='.')))
fdata <- new("AnnotatedDataFrame", data=data.frame(matrix(z$dnew,nrow=nsel,byrow=TRUE),matrix(z$anew,nrow=nsel,byrow=TRUE),matrix(z$lnew,nrow=nsel,byrow=TRUE)),varMetadata=metadata)
sampleNames(fdata) <- paste("Gene",1:nrow(fdata))
varLabels(fdata) <- c(paste('d',1:length(groupsnewr),sep='.'),paste('alpha',1:length(groupsnewr),sep='.'),paste('mean',1:length(groupsnewr),sep='.'))

experimentData <- new("MIAME", title = "Dataset simulated with simnewsamples", abstract = "This dataset contains expression data simulated from the posterior predictive distribution of a GaGa or MiGaGa. The expression data can be accessed via exprs(object), and the parameter values used to generate the data through fData(object)")

x <- data.frame(matrix(z$xnew,nrow=nsel,byrow=TRUE))
names(x) <- paste('Array',1:nrow(pheno)); rownames(x) <- paste('Gene',1:nrow(fdata))
x <- new("ExpressionSet", phenoData=pheno, featureData=fdata, exprs = x, experimentData = experimentData)
return(x)

}
