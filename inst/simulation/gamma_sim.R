library(Matrix)
library(mvtnorm)
library(ncvreg)
source("../sample.splitting.dim.reduction.R")
source("../PS5_mediation.R")
gamma.sim <- function(
    N, 
    p,
    cov, 
    num.ts,
    replication = 100
){
  alpha.max <- ifelse(num.ts<=10, 0.4, 0.2)
  Alpha.effect.setting <- seq(0, alpha.max, length.out = 21)
  Alpha.effect.setting <- Alpha.effect.setting[-1]
  power.res <- data.frame(matrix(NA, ncol = 2, nrow = length(Alpha.effect.setting), 
                                 dimnames=list(NULL, c("PS5.L1", "PS5.L2"))))
  for(kk in 1:length(Alpha.effect.setting)){
    # power setting
    Beta.effect <- 0.5
    Alpha.c <- rep(0.0, p)								                            # C -> M (the association between covariate and mediators)
    Beta.x <- 0.3 									                                  # X -> Y (the association between exposure and outcome [conditional on mediators])
    Beta.c <- 0.1                                                     # C -> Y (the association between covariate and outcome)
    Alpha.x <- rep(0, p)
    Alpha.x[1:num.ts] <- Alpha.effect.setting[kk]
    Beta.m <- c(rep(Beta.effect, 50), rep(0, p-50))
    if(cov == 0.9){
      # create moderate correlation structure
      exponent <- abs(matrix(1:p - 1, nrow = p, ncol = p, byrow = TRUE) - (1:p - 1))
      v <- 0^exponent
      # ultra-high correlation among pair of mediators (e.g. M1 and M2, M3 and M4)
      for(i in 1:floor(num.ts/2)){
        v[(2*i-1),2*i] <- v[2*i,(2*i-1)] <- cov
      }
    }else{
      exponent <- abs(matrix(1:p - 1, nrow = p, ncol = p, byrow = TRUE) - (1:p - 1))
      v <- cov^exponent                                                 # covariance matrix of residual (AR1 model)
    }
    
    # start
    All.PVAL <- data.frame(matrix(NA, ncol = 2, nrow = replication, 
                                  dimnames=list(NULL, c("PS5.L1", "PS5.L2"))))
    for(rep.time in 1:replication){
      print(rep.time)
      set.seed(rep.time*1357)
      epsilon <- rmvnorm(n=N, mean=rep(0, p), as.matrix(v))
      exposure <- rnorm(N, mean = 0, sd = 1)
      conf <- c(rep(0, ceiling(N/6)), rep(1, floor(N/3)), rep(0, floor(N/3)), rep(1, ceiling(N/6)))
      conf2<-runif(N, 0, 1)
      mediator <- 0.3 + (conf%o%Alpha.c+exposure%o%Alpha.x+epsilon)
      colnames(mediator) <- paste0("M", c(1:p))
      error <- rnorm(N, 0, 1)
      diz <- 0.2 + conf*Beta.c + exposure*Beta.x + mediator%*%Beta.m + error  # no interaction
      conf.mat <- matrix(cbind(conf,conf2), ncol = 2)
      exposure.mat <- matrix(exposure, ncol = 1)
      # L1
      PS5.L1 <- PS5.mediation(M = mediator, 
                              X = exposure.mat, 
                              Y = diz, 
                              C = conf.mat,
                              n.draw = 10000,
                              dim.reduction = T,
                              ga = 1)
      
      # L2
      PS5.L2 <- PS5.mediation(M = mediator, 
                              X = exposure.mat, 
                              Y = diz, 
                              C = conf.mat,
                              n.draw = 10000,
                              dim.reduction = T,
                              ga = 2)
      
      ######################################################################
      ############################### power ################################ 
      ######################################################################
      All.PVAL[rep.time,] <- c(PS5.L1$global.test, 
                               PS5.L2$global.test)
    }
    # output 
    power.res[kk,] <- apply(All.PVAL<0.05, 2, function(x) mean(x, na.rm=T))
  }
  return(power.res)
}

## simulation
cov.setting <- 0       # c(0,0.5,0.9)
num.ts.setting <- 1    # c(5,10,30)
N.setting <- 500
p.setting <- 1000
print(paste0('num.ts = ', num.ts.setting, '; cor = ', cov.setting))
gamma.result <- gamma.sim(N = N.setting,
                          p = p.setting,
                          cov = cov.setting,
                          num.ts = num.ts.setting,
                          replication = 100)
save(gamma.result, file = paste0("gamma_N", N.setting, "p", p.setting,"_ts", num.ts.setting, ".RData"))


