simnewsamples.gagafit <- function(gg.fit,groupsnew,sel,x,groups) {
# Simulate parameters from the posterior of a GaGa or MiGaGa model, and new observations from the posterior predictive

gapprox <- TRUE
if (is(x, "exprSet") | is(x,"ExpressionSet")) {
  if (is.character(groups) && length(groups)==1) { groups <- as.factor(pData(x)[, groups]) }
  x <- exprs(x)
} else if (!is(x,"data.frame") & !is(x,"matrix")) { stop("x must be an exprSet, data.frame or matrix") } 

patterns <- gg.fit$patterns
v <- gg.fit$pp

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
par <- getpar(gg.fit)
alpha0 <- as.double(par$a0); nu <- as.double(par$nu); balpha <- as.double(par$balpha)
nualpha <- as.double(par$nualpha)
nclust <- as.integer(gg.fit$nclust); rho <- as.double(par$probclus)
sumx <- double(nrow(x)*sum(ngrouppat)); nobsx <- double(sum(ngrouppat))
if (gg.fit$equalcv) {
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

z <- .C("compute_sumxC",sumx=t(sumx),prodx=t(prodx),nobsx=nobsx,as.integer(gg.fit$equalcv),nsel,sel,as.integer(sum(ngrouppat)),ncol(x),as.double(t(x)),groupsr,K,as.integer(nrow(patterns)),as.integer(t(patterns)),ngrouppat,as.integer(1))
sumx <- matrix(z$sumx,nrow=nrow(x),byrow=TRUE); prodx <- matrix(z$prodx,nrow=nrow(x),byrow=TRUE)
nobsx <- z$nobsx

z <- .C("simnewsamples_ggC",xnew=xnew,dnew=dnew,anew=anew,lnew=lnew,nsamples,groupsnewr,K,nsel,sel,alpha0,nu,balpha,nualpha,as.integer(gg.fit$equalcv),nclust,rho,as.double(t(v)),as.integer(nrow(patterns)),as.integer(t(patterns)),ngrouppat,t(sumx),t(prodx),nobsx,gapprox)

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
