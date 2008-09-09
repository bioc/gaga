checkfit.gagafit <- function(gg.fit,x,groups,type='data',logexpr=FALSE,xlab,ylab,main,lty,lwd,...) {
# Plots to check fit of GaGa and MiGaGa models.

if (is(x, "exprSet") | is(x, "ExpressionSet")) {
  if (is.character(groups) && length(groups)==1) { groups <- as.factor(pData(x)[, groups]) }
  x <- exprs(x)
} else if (!is(x,"data.frame") & !is(x,"matrix")) { stop("x must be an ExpressionSet, exprSet, data.frame or matrix") }

if (ncol(x) != length(groups)) stop('The length of argument groups does not match with the number of columns in x')
if ((type!='data') && (type!='shape') && (type!='mean') && (type!='shapemean')) stop('The argument type is not valid')

xpred <- simnewsamples(gg.fit=gg.fit,groupsnew=groups,x=x,groups=groups)
if (type=='data') {
  if (logexpr) {
    xnewpdf <- density(log2(unlist(exprs(xpred)))); if (is.list(x)) xpdf <- density(log2(unlist(x))) else xpdf <- density(log2(x))
  } else {
    xnewpdf <- density(unlist(exprs(xpred))); if (is.list(x)) xpdf <- density(unlist(x)) else xpdf <- density(x)
  }
  if (missing(xlab)) xlab <- 'Expression values'; if (missing(ylab)) ylab <- 'Density'; if (missing(main)) main <- ''
  plot(xpdf,type='l',xlab=xlab,ylab=ylab,main=main,...); lines(xnewpdf,lty=2,lwd=2); legend(max(xpdf$x),max(xpdf$y),c('Observed data','Posterior predictive'),lty=1:2,lwd=1:2,xjust=1,yjust=1) 
} else if (type=='shape') {
  colsel <- (1+ncol(fData(xpred))/3):(2*ncol(fData(xpred))/3)
  aest <- rowMeans(x)^2/apply(x,1,'var'); apdf <- density(unlist(fData(xpred)[,colsel]))
  if (missing(xlab)) xlab <- 'alpha parameters (shape)'; if (missing(ylab)) ylab <- 'Density'; if (missing(main)) main <- ''
  plot(apdf,xlab=xlab,ylab=ylab,main=main,lty=1,...); lines(density(aest),lty=2)
  legend(quantile(aest,probs=.95),max(apdf$y),c('Model-based','Moments estimate'),lty=1:2,xjust=0,yjust=1)
} else if (type=='mean') {
  colsel <- (1+2*ncol(fData(xpred))/3):ncol(fData(xpred))
  lest <- rowMeans(x); lpdf <- density(unlist(fData(xpred)[,colsel]))
  if (missing(xlab)) xlab <- 'mean parameters'; if (missing(ylab)) ylab <- 'Density'; if (missing(main)) main <- ''
  plot(lpdf,xlab=xlab,ylab=ylab,main=main,lty=1,...); lines(density(lest),lty=2)
  legend(max(lest),max(lpdf$y),c('Model-based','Moments estimate'),lty=1:2,xjust=1,yjust=1)
} else if (type=='shapemean') {
  colsel1 <- (1+ncol(fData(xpred))/3):(2*ncol(fData(xpred))/3)
  colsel2 <- (1+2*ncol(fData(xpred))/3):ncol(fData(xpred))
  aest <- rowMeans(x)^2/apply(x,1,'var'); lest <- rowMeans(x)
  if (missing(xlab)) xlab <- 'mean parameters'; if (missing(ylab)) ylab <- 'alpha parameters (shape)'; if (missing(main)) main <- ''
  plot(unlist(fData(xpred)[,colsel2]),unlist(fData(xpred)[,colsel1]),xlab=xlab,ylab=ylab,main=main,...); points(lest,aest,col=2)
}

}
