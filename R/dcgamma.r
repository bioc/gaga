dcgamma <- function(x,a,b,c,d,r,s,newton=TRUE) {
# Returns density of a conjugate gamma shape distribution evaluated at x
# Input:
# - x: vector with points to evaluate the density at
# - a,b,c,d,r,s,p: parameters
#  - newton: newton==TRUE tries a Newton step to improve Gamma approximation
# Output: 
# - y: density of a conjugate gamma shape distribution with parameters a,b,c evaluated at x

if ((a<0) | (b<0) | (d<0) | (r<0) | (s<0)) stop('Parameters a,b,d,r,s must be >=0')
if (a==0) { 
  if (b-d<=0) stop('Non-valid parameters. b must be > than (p+1)*d')
  if (c<=0) stop('Non-valid parameters. c must be > 0')
} else {
  if (b+.5*a-.5<=0) stop('Non-valid parameters. b must be > .5-.5*a')
  if (c+a*log(s/a)<=0) stop('Non-valid parameters. c must be > a*log(a/s)')
}

x <- as.double(x); n <- as.integer(length(x)); a <- as.double(a); b <- as.double(b); c <- as.double(c); d <- as.double(d); r <- as.double(r); s <- as.double(s); y <- double(n); newton <- as.integer(newton)
normk <- as.double(-1)            #indicates the C routine to calculate the normalization constant
z <- .C("dcgammaC",y=y,normk,x,n,a,b,c,d,r,s,newton)
return(z$y)
}
