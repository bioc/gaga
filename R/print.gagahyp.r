print.gagahyp <- function(patterns) {
  groups <- paste("Group",1:ncol(patterns))
  for (i in 1:nrow(patterns)) {
    cnt <- length(table(patterns[i,]))
    cat("  Pattern ",i-1,": ",sep="")
    for (j in 0:(cnt-1)) {
      if (j<(cnt-1)) { eq <- c(rep(" =",sum(patterns[i,]==j)-1),"!=") } else { eq <- c(rep("=",sum(patterns[i,]==j)-1),"\n") }
      cat(paste(groups[patterns[i,]==j],eq))
    }
  }
}
