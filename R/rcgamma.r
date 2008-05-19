rcgamma <- function(n,a,b,c,d,r,s,newton=TRUE) {
# Random draws from a conjugate gamma shape distribution by approximating it with a Gamma

if ((sum(a<0)>0) | (b<0) | (d<0) | (r<0) | (sum(s<0)>0)) stop('Parameters a,b,d,r,s must be >=0')
if ((sum(a)+sum(s))==0) { 
  if (b<=0) stop('Non-valid parameters. b must be > 0')
  if (c<=0) stop('Non-valid parameters. c must be > 0')
} else {
  if (sum(b+.5*a-1.5)+1<=0) stop('Non-valid parameters. sum(b+.5*a-1.5)+1 must be > 0')
  if (sum(c+a*log(s/a))<=0) stop('Non-valid parameters. sum(c+a*log(s/a)) must be > 0')
}
if (length(a)!=length(s)) stop('Arguments a and s must be vectors of the same length')

n <- as.integer(n); a <- as.double(a); b <- as.double(b); c <- as.double(c); d <- as.double(d); r <- as.double(r); s <- as.double(s)
newton <- as.integer(newton)

x <- double(n)
z <- .C("rcgammaC",x=x,n,a,b,c,d,r,s,as.integer(length(a)),newton)
return(z$x)

}
