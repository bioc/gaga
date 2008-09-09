mcgamma <- function(a,b,c,d,r,s,newton=TRUE) {
# Moments for a conjugate gamma shape distribution */

if ((sum(a<0)>0) | (b<0) | (d<0) | (r<0) | (sum(s<0)>0)) stop('Parameters a,b,d,r,s must be >=0')
if ((sum(a)+sum(s))==0) { 
  if (b<=0) stop('Non-valid parameters. b must be > 0')
  if (c<=0) stop('Non-valid parameters. c must be > 0')
} else {
  if (b+sum(.5*a-1.5)+1<=0) stop('Non-valid parameters. b+sum(.5*a-1.5)+1 must be > 0')
  if (c+sum(a*log(s/a))<=0) stop('Non-valid parameters. c+sum(a*log(s/a)) must be > 0')
}
if (length(a)!=length(s)) stop('Arguments a and s must be vectors of the same length')

a <- as.double(a); b <- as.double(b); c <- as.double(c); d <- as.double(d); r <- as.double(r); s <- as.double(s)
newton <- as.integer(newton)
normk <- m <- as.double(-1); v <- double(1)
z <- .C("mcgammaC",normk=normk,m=m,v=v,a,b,c,d,r,s,as.integer(length(a)),newton)
return(list(m=z$m,v=z$v,normk=z$normk))

}
