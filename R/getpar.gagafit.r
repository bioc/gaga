getpar.gagafit <- function(fit) {
#returns parameter estimates from a 'gagafit' object in a named list
if (sum(is.na(fit$parest)>0)) stop('Parameter estimates not available. Use function parest first.')
nclust <- fit$nclust
return(list(a0=fit$parest[1:nclust],nu=fit$parest[(nclust+1):(2*nclust)],balpha=fit$parest[2*nclust+1],nualpha=fit$parest[2*nclust+2],probclus=fit$parest[(2*nclust+3):(3*nclust+2)],probpat=fit$parest[-1:-(3*nclust+2)]))
}
