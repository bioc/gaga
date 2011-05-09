buildPatterns <- function(groups) {
  getPatterns <- function(x) {
    tmpList <- vector('list',length=nrow(x))
    for (i in 1:nrow(x)) {
      mymax <- max(x[i,])
      tmp <- matrix(ncol=ncol(x)+1,nrow=mymax+2)
      for (j in 1:nrow(tmp)) tmp[j,] <- c(x[i,],j-1)
      tmpList[[i]] <- tmp
    }
    do.call('rbind',tmpList)
  }
  ans <- matrix(0,ncol=1,nrow=1)
  for (i in 2:length(unique(groups))) {
#    cat(paste('Pattern',i,'\n'))
    ans <- getPatterns(ans)
  }
  colnames(ans) <- as.character(unique(groups))
  return(ans)
}
