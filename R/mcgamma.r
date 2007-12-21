mcgamma <- function(a,b,c,d,r,s,newton=TRUE) {
# Computes moments for a conjugate gamma shape distribution */
# Input:
#  - a,b,c,d,r,s: parameters
#  - newton: newton==TRUE tries a Newton step to improve Gamma approximation
# Output: list with the following elements
#  - m: mean 
#  - v: variance 
#  - normk: normalization constant

if ((a<0) | (b<0) | (d<0) | (r<0) | (s<0)) stop('Parameters a,b,d,r,s must be >=0')
if (a==0) { 
  if (b-d<=0) stop('Non-valid parameters. b must be > than d')
  if (c<=0) stop('Non-valid parameters. c must be > 0')
} else {
  if (b+.5*a-.5<=0) stop('Non-valid parameters. b must be > .5-.5*a')
  if (c+a*log(s/a)<=0) stop('Non-valid parameters. c must be > a*log(a/s)')
}

a <- as.double(a); b <- as.double(b); c <- as.double(c); d <- as.double(d); r <- as.double(r); s <- as.double(s); newton <- as.integer(newton)
normk <- m <- as.double(-1); v <- double(1)
z <- .C("mcgammaC",normk=normk,m=m,v=v,a,b,c,d,r,s,newton)
return(list(m=z$m,v=z$v,normk=z$normk))

}
