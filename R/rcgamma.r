rcgamma <- function(n,a,b,c,d,r,s,newton=TRUE) {
# Generates random draws from a conjugate gamma shape distribution by approximating it with a Gamma */
# Input:
# - n: number of random draws to generate
# - a,b,c,e,r,s: parameters
#  - newton: newton==TRUE tries a Newton step to improve Gamma approximation
# Output:
# - x: vector of length n with the random draws

if ((a<0) | (b<0) | (d<0) | (r<0) | (s<0)) stop('Parameters a,b,d,r,s must be >=0')
if (a==0) { 
  if (b-d<=0) stop('Non-valid parameters. b must be > than d')
  if (c<=0) stop('Non-valid parameters. c must be > 0')
} else {
  if (b+.5*a-.5<=0) stop('Non-valid parameters. b must be > .5-.5*a')
  if (c+a*log(s/a)<=0) stop('Non-valid parameters. c must be > a*log(a/s)')
}

n <- as.integer(n); a <- as.double(a); b <- as.double(b); c <- as.double(c); d <- as.double(d); r <- as.double(r); s <- as.double(s); newton <- as.integer(newton)

x <- double(n)
z <- .C("rcgammaC",x=x,n,a,b,c,d,r,s,newton)
return(z$x)

}
