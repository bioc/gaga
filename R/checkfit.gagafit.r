checkfit.gagafit <- function(gg.fit,x,groups,type='data',...) {
# Plots to check fit of TripleG model. Compares posterior of parameters with method of moments estimates or observed data with posterior predictive. To compare with prior use function sim.gg.
# - x: data used to fit the model
# - groups: vector of length ncol(x) indicating what group does each column in x correspond to
# - gg.fit: model fit with parameter estimates as returned by parest.gagafit
# - type: 'data' for marginal density of the data; 'shape' for shape parameter; 'mean' for mean parameter; 'shapemean' for joint of shape and mean parameters
# - ...: other arguments to be passed to the plot function

if (ncol(x) != length(groups)) stop('The length of argument groups does not match with the number of columns in x')
if ((type!='data') && (type!='shape') && (type!='mean') && (type!='shapemean')) stop('The argument type is not valid')

xpred <- simnewsamples(gg.fit=gg.fit,groupsnew=groups,x=x,groups=groups)
if (type=='data') {
  xnewpdf <- density(xpred$xnew); xpdf <- density(x)
  plot(xpdf,type='l',...); lines(xnewpdf,lty=2,lwd=2); legend(max(xpdf$x),max(xpdf$y),c('Observed data','Posterior predictive'),lty=1:2,lwd=1:2,xjust=1,yjust=1) 
} else if (type=='shape') {
  aest <- rowMeans(x)^2/apply(x,1,'var'); apdf <- density(xpred$anew)
  plot(apdf,lty=1,...); lines(density(aest),lty=2)
  legend(quantile(aest,probs=.95),max(apdf$y),c('Model-based','Moments estimate'),lty=1:2,xjust=0,yjust=1)
} else if (type=='mean') {
  lest <- rowMeans(x); lpdf <- density(xpred$lnew)
  plot(lpdf,lty=1,...); lines(density(lest),lty=2)
  legend(max(lest),max(lpdf$y),c('Model-based','Moments estimate'),lty=1:2,xjust=1,yjust=1)
} else if (type=='shapemean') {
  aest <- rowMeans(x)^2/apply(x,1,'var'); lest <- rowMeans(x)
  plot(xpred$lnew,xpred$anew,...); points(lest,aest,col=2)
}

}
