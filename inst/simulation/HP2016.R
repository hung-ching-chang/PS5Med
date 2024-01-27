# HP2016
HP2016 <-function(
    M = mediators,	      	# n-by-p matrix
    X = exposure,		      	# n-by-1 matrix
    Y = outcome,		      	# n-by-1 matrix
    C = conf,		          	# n-by-q matrix, not including intercept
    n.draw = 1000,	      	# no. of Monte-Carlo resampling
    adaptive.ratio = 0.8    # dimension reduction
){
  set.seed(1004)
  n <- dim(M)[1]
  p <- dim(M)[2]

  ### transformation
  fit.m <- lm(M~C+X)
  Sigma <- cov(fit.m$residual) # covariance matrix of error term
  svds <- tryCatch(svd(Sigma),
                   error=function(e) NULL)
  if(is.null(svds)){return(list(HP2016 = NA,
                                global.me = NA))}
  if (adaptive.ratio != 1){
    n.factor<-sum(cumsum(svds$d)/sum(svds$d)<adaptive.ratio)+1
    U<-svds$u[,1:n.factor]
    p<-n.factor
  } else {
    U<-svds$u
  }
  M.star <- M%*%U

  ### estimation
  fit.m2 <- lm(M.star~C+X)
  if(p>1){
    cov1.star <- diag(diag(vcov(fit.m2))[(1:p)*(2+ncol(C))])  # (2+ncol(C)) to select cov of alpha star
  }else{
    cov1.star <- diag(vcov(fit.m2))[(1:p)*(2+ncol(C))]
  }
  beta.star <- NULL   # beta star estimation
  mse <- NULL   # mse for sigma^2 estimation
  for (j in 1:p){
    M.star.j <- M.star[,j]
    fit.y2Mj <- lm(Y~C+X+M.star.j)
    beta.star <- c(beta.star, fit.y2Mj$coef[c("M.star.j")])
    mse <- c(mse, summary(fit.y2Mj)$sigma^2)
  }
  if(p>1){
    alpha.star <- fit.m2$coef["X",]   # alpha star estimation
  }else{
    alpha.star <- fit.m2$coef["X"]   # alpha star estimation
  }
  me.star <- beta.star*alpha.star # mediation effect = alpha * beta

  ### covariance
  Z <- cbind(1, X, C, M.star)
  info <- t(Z)%*%Z
  cov2.star <- min(mse) * solve(info)[2+ncol(C)+1:p, 2+ncol(C)+1:p]
  bcov <- bdiag(cov1.star, cov2.star)  # block diagonal matrix

  ### L2 norm method (Huang and Pan 2016)
  theta.mc <- rmvnorm(n.draw,
                      mean = c(alpha.star, beta.star),
                      sigma = as.matrix(bcov))
  alpha.mc <- theta.mc[,1:p]
  beta.mc <- theta.mc[,p+1:p]
  me.mc <- as.matrix(alpha.mc*beta.mc)   # alpha * beta
  me.mc.c <- t(t(me.mc) - apply(me.mc, 2, mean))
  T.mc <- apply(me.mc.c^2, 1, sum)	# T statistic
  pval.L2 <- mean(T.mc > sum(me.star^2))

  return(list(HP2016 = pval.L2,
              global.me = sum(me.star)))
}


