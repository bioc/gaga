getpar.gagafit <- function(gg.fit) {
#returns parameter estimates from a 'gagafit' object in a named list
if (sum(is.na(gg.fit$parest)>0)) stop('Parameter estimates not available. Use function parest first.')
nclust <- gg.fit$nclust
return(list(a0=gg.fit$parest[1:nclust],nu=gg.fit$parest[(nclust+1):(2*nclust)],balpha=gg.fit$parest[2*nclust+1],nualpha=gg.fit$parest[2*nclust+2],probclus=gg.fit$parest[(2*nclust+3):(3*nclust+2)],probpat=gg.fit$parest[-1:-(3*nclust+2)]))
}
