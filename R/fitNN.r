fitNN <- function(x, groups, patterns, B=20, trace=TRUE) {
  require(EBarrays)
  if (is(x, "exprSet") | is(x, "ExpressionSet")) {
    if (is.character(groups) && length(groups) == 1) {
      groups <- as.factor(pData(x)[, groups])
    }
    x <- exprs(x)
  }
  else if (!is(x, "data.frame") & !is(x, "matrix")) {
      stop("x must be an ExpressionSet, exprSet, data.frame or matrix")
  }
  #Patterns
  K <- length(unique.default(groups))
  if (missing(patterns)) {
    patterns <- rbind(rep(0, K), 0:(K - 1))
    colnames(patterns) <- unique.default(groups)
  }
  groupsr <- groups2int(groups, patterns) + 1
  if (length(table(groupsr)) != ncol(patterns)) stop("patterns must have the same number of columns as the number of distinct groups")
  if (ncol(patterns) != K) stop("patterns must have number of columns equal to the number of distinct elements in groups")
  if (sum(is.na(patterns)) > 0) stop("patterns cannot have any NA values")
  if (sum(is.nan(patterns)) > 0) stop("patterns cannot have any NaN values")
  if (sum(is.infinite(patterns)) > 0) stop("patterns cannot have any Inf values")
  for (i in 1:nrow(patterns)) { patterns[i, ] <- as.integer(as.integer(as.factor(patterns[i,])) - 1) }
  class(patterns) <- "gagahyp"
  hypotheses <- makeEBarraysHyp(patterns=patterns, groups=groupsr)
  #Fit model
  expx <- exp(x)
  verbose <- getOption("verbose")
  options(verbose=trace)
  family <- eb.createFamilyLNNMV()
  nn.fit <- emfit(data=expx, family=family, hypotheses=hypotheses, groupid=groupsr, num.iter=B)
  options(verbose=verbose)
  priorest <- sigmaPriorEst(x=x,groupid=groupsr,model='NN')
  pp <- postprob(fit=nn.fit, data=expx, groupid=groupsr)$pattern
  #Return output
  parest <- c(mu0=nn.fit@thetaEst[1,'theta1'],tau02=exp(nn.fit@thetaEst[1,'theta2']),priorest['v0'],priorest['sigma02'],probclus=1,probpat=nn.fit@probEst)
  ans <- list(parest=parest, patterns=patterns, pp=pp, nn.fit=nn.fit)
  class(ans) <- 'nnfit'
  return(ans)
}



sigmaPriorEst <- function(x, groupid, model='LNN') {
#Estimate sigma_j^2 ~ IG(.5*v0,.5*v0*sigma0^2) prior LNN & NN models via Method of Moments (as in package EBarrays)
  if (is(x, "ExpressionSet")) x <- exprs(x)
  if (model=='LNN') { x <- log(x) } else if (model!='NN') { stop('model not recognized') }
  n <- sum(groupid != 0)
  groups <- unique(groupid)
  groups <- groups[groups != 0]
  ngroup <- length(groups)
  vars <- 0
  for (i in 1:ngroup) {
    temp <- x[, groupid == groups[i]]
    ncoltemp <- sum(groupid == groups[i])
    if (ncoltemp == 1) { vars <- vars } else { vars <- vars + (ncoltemp - 1) * apply(temp, 1, FUN = var) }
  }
  svars <- vars/(n - ngroup)
  v0 <- mean(svars)^2/((length(svars) - 1)/length(svars) * var(svars)) * 2 + 4
  sigma02 <- mean(svars) * (v0 - 2)/v0
  return(c(v0=v0,sigma02=sigma02))
}

#Format patterns & groups information into ebarraysPatterns object. Used to call EBarrays routines
makeEBarraysHyp <- function(patterns, groups) {
hypotheses <- character(nrow(patterns))
for (i in 1:nrow(patterns)) {
  tmp <- integer(length(groups))
  for (j in unique(patterns[i,])) tmp[groups %in% colnames(patterns)[patterns[i,]==j]] <- j+1
  hypotheses[i] <- paste(tmp,collapse=' ')
}
#  ngroup <- table(groups)
#  for (i in 1:length(hypotheses)) hypotheses[i] <- paste(rep(patterns[i,]+1,ngroup),collapse=' ')
  ebPatterns(hypotheses)
}
