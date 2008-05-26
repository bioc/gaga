print.gagahyp <- function(x,probpat=NA,...) {
  groups <- colnames(x)
  for (i in 1:nrow(x)) {
    cnt <- length(table(x[i,]))
    cat("  Pattern ",i-1," ",sep="")
    if (sum(is.na(probpat))==0) {
      if (is.integer(probpat)) {
        cat("(",probpat[i]," genes): ",sep='') 
      } else {
        cat("(",round(100*probpat[i],1),"% genes): ",sep='')
      }
    } else { cat(": ") }
    for (j in 0:(cnt-1)) {
      if (j<(cnt-1)) { eq <- c(rep(" =",sum(x[i,]==j)-1),"!=") } else { eq <- c(rep("=",sum(x[i,]==j)-1),"\n") }
      cat(paste(groups[x[i,]==j],eq))
    }
  }
}
