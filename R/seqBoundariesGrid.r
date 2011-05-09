seqBoundariesGrid <- function(b0,b1,forwsim,samplingCost,powmin=0,f='linear',ineq='less') {
# Expected utility on a grid of parametric sequential boundaries
 
if ((!is.vector(b0)) | (!is.vector(b1))) stop('b0 and b1 must be vectors')
if (missing(forwsim)) stop('forwsim must be specified')
if (missing(samplingCost)) stop('samplingCost must be specified')
if ((f!='linear') & (f!='invsqrt')) stop('f must be linear or invsqrt')
if ((ineq!='greater') & (ineq!='less')) stop('ineq must be greater or less')

forwsim$summary[is.na(forwsim$summary)] <- -9999  #indicate NAs with -9999
b0opt <- b1opt <- uopt <- fdropt <- fnropt <- poweropt <- jopt <- double(1)
nb0 <- as.integer(length(b0)); nb1 <- as.integer(length(b1)); nsim <- as.integer(length(unique(forwsim$simid)))
ave <- TRUE
if (ave) { ustop <- jstop <- double(nb0*nb1) } else { ustop <- jstop <- double(nb0*nb1*nsim) }
fdrstop <- fnrstop <- powerstop <- double(nb0*nb1)
minj <- as.double(min(forwsim$time))
f <- ifelse(f=='linear',1,2); ineq <- ifelse(ineq=='greater',1,-1)
J <- as.integer(max(forwsim$time)-minj)

b0grid <- b1grid <- double(nb0*nb1)
optlast <- as.integer(1); fdrmax <- as.double(1); fnrmax <- double(1) #currently unused variables (passed to C)

ans <- .C("euCgrid",b0opt=b0opt,b1opt=b1opt,uopt=uopt,fdropt=fdropt,fnropt=fnropt,poweropt=poweropt,jopt=jopt,b0grid=b0grid,b1grid=b1grid,ustop=ustop,fdrstop=fdrstop,fnrstop=fnrstop,powerstop=powerstop,jstop=jstop,as.double(samplingCost),nb0,nb1,as.double(b0),as.double(b1),as.integer(ave),nsim,as.integer(nrow(forwsim)),as.integer(forwsim$simid),as.double(forwsim$time),minj,as.double(forwsim$u),as.double(forwsim$fdr),as.double(forwsim$fnr),as.double(forwsim$power),as.double(forwsim$summary),as.double(fdrmax),as.double(fnrmax),as.double(powmin),as.integer(f),as.double(ineq),J,optlast)

opt <- c(b0=ans$b0opt,b1=ans$b1opt,u=ans$uopt,fdr=ans$fdropt,fnr=ans$fnropt,power=ans$poweropt,time=ans$jopt)
grid <- data.frame(b0=ans$b0grid,b1=ans$b1grid,u=ans$ustop,fdr=ans$fdrstop,fnr=ans$fnrstop,power=ans$powerstop,time=ans$jstop)

return(list(opt=opt,grid=grid))
}
