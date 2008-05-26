print.gagaclus <- function(x,...) {
  print(x$patterns,table(factor(x$d,levels=0:nrow(x$patterns))))
}
